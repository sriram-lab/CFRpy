"""
This module provides methods for the Python implementation of the constrain flux regulation (CFR) algorithm.

"""

__licence__ = "GPL GNU"
__docformat__ = "reStructuredText"

# Dependencies
import os, sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import gurobipy as gp
from gurobipy import GRB, GurobiError
import scipy.io as sio
import scipy.sparse as sp
from itertools import chain
from warnings import warn
from cobra.core import Solution
from cobra.io.mat import from_mat_struct
from cobra.util.array import create_stoichiometric_matrix
from cobra.manipulation.delete import knock_out_model_genes
from optlang.glpk_interface import Variable, Constraint, Objective, Model

# Methods
def load_matlab_model_rxnGeneMat(infile_path, variable_name=None, inf=np.inf):
    """
    Load a cobra model from a .mat file, identical to
    cobra.io.load_matlab_model, but also attaches model.rxnGeneMat: a
    reactions x genes scipy.sparse matrix (row/col order matching
    model.reactions/model.genes at load time), or None if the .mat file
    has no rxnGeneMat field.

    'rxnGeneMat' is not a name any cobrapy method reads or writes, so this
    attribute is purely additive -- it rides along on the Model object
    without changing the behavior of any other cobrapy function:
      - Model has no __slots__, so arbitrary attributes are allowed.
      - model.copy() deep-copies whatever is on the object, including this.
      - save_matlab_model always rebuilds its own rxnGeneMat fresh from
        gene_reaction_rule on export -- it never reads model.rxnGeneMat,
        so round-tripping through save/load is unaffected either way.

    Caveat: the row/column alignment only holds at load time. If
    reactions or genes are later added, removed, or reordered, this
    attribute becomes stale and is not kept in sync automatically.
    """
    data = sio.loadmat(infile_path)
    meta_vars = {"__globals__", "__header__", "__version__"}

    if variable_name is not None:
        model = from_mat_struct(data[variable_name], model_id=variable_name, inf=inf)
        m = data[variable_name]
    else:
        possible_names = sorted(k for k in data if k not in meta_vars)
        model, m = None, None
        for name in possible_names:
            try:
                model = from_mat_struct(data[name], model_id=name, inf=inf)
                m = data[name]
                break
            except ValueError as e:
                print(f"Some problem with the model, causing error {e}")
        if model is None:
            raise IOError(f"No COBRA model found at {infile_path}.")

    if "rxnGeneMat" in (m.dtype.names or ()):
        rgm = m["rxnGeneMat"][0, 0]
        if not sp.issparse(rgm):
            rgm = sp.csr_matrix(rgm)
        model.rxnGeneMat = rgm.tocsr()
    else:
        model.rxnGeneMat = None

    return model


def geneKO_from_rxnGeneMat(cobra_model, gene_list):
    """
    Mirrors cobrapy's knock_out_model_genes, but determines the candidate
    reaction set (rxn_set) from gene_rxn_map -- built from a model's
    rxnGeneMat -- instead of the model's own live gene.reactions.

    Genes are still knocked out via gene.knock_out(), so the functional
    evaluation itself is unchanged; only *discovery* of which reactions
    to check is sourced externally. Call this inside `with cobra_model:`
    to keep the knockouts (and their bound changes) transient.
    """
    rgm = getattr(cobra_model, 'rxnGeneMat', None)
    rgm = rgm.tocsc()
    gene_ids = [g.id for g in cobra_model.genes]
    rxn_ids = [r.id for r in cobra_model.reactions]
    gene_rxn_map = {
        gene_ids[j]: [rxn_ids[i] for i in rgm.indices[rgm.indptr[j]:rgm.indptr[j + 1]]]
        for j in range(rgm.shape[1])
    }
    rxn_set = set()
    for gene in cobra_model.genes.get_by_any(gene_list):
        gene.knock_out()
        for rxn_id in gene_rxn_map.get(gene.id, []):
            try:
                rxn_set.add(cobra_model.reactions.get_by_id(rxn_id))
            except KeyError:
                continue  # rxn_id in the map doesn't exist in this cobra model
    return [rxn for rxn in rxn_set if not rxn.functional]


def opt_variable(name, lb, ub, vtype): 
    """
    Define an optlang variable.

    Parameters
    ----------
    name : str
        A string specifying the variable name.
    lb : numeric
        The numeric value assigned to the lower bound.
    ub : numeric
        The numeric value assigned to the upper bound.
    vtype : single-character string
        Determines the variable type (B = binary, C = continuous).

    Returns
    -------
    Variable
        A variable object with defined type and bounds. 

    Examples
    --------
    Usage cases of the `opt_variable` function. 
    
    Case 1 (binary variable): 
    >>> name, lb, ub, vtype = 'x', 0, 1, 'B'
    >>> x = opt_variable(name, lb, ub, vtype)
    >>> print(x)
    0 <= x <= 1

    Case 2 (continuous variable): 
    >>> name, lb, ub, vtype = 'x', -1.5, 1.5, 'C'
    >>> x = opt_variable(name, lb, ub, vtype)
    >>> print(x)
    -1.5 <= x <= 1.5

    """

    # Define variable
    if vtype=='B': 
        x = Variable(name, lb=lb, ub=ub, type='binary')
    else: 
        x = Variable(name, lb=lb, ub=ub, type='continuous')

    # Return output
    return x

def opt_constraint(x, sense, b): 
    """
    Define an optlang constraint.
    
    Parameters
    ----------
    x : Variable
        An optlang variable to constrain.
    sense : single-character string (choose among '<', '>', '=')
        Specifies the comparative relationship between the left-hand and right-hand side of the constraint equation.
    b : numeric
        The numeric value assigned to the right-hand side (i.e., bound) of the constraint expression.

    Returns
    -------
    Constraint
        A constraint object with defined type and bounds. 

    Examples
    --------
    Usage cases of the `opt_constraint` function. 
    
    Case 1 (greater than): 
    >>> x, sense, b = x, '>', -10
    >>> c = opt_constraint(x, sense, b)
    >>> print(c)
    -10 <= x

    Case 2 (less than): 
    >>> x, sense, b = x, '<', 10
    >>> c = opt_constraint(x, sense, b)
    >>> print(c)
    x <= 10

    Case 3 (equal): 
    >>> x, sense, b = x, '=', 10
    >>> c = opt_constraint(x, sense, b)
    >>> print(c)
    10 <= x <= 10

    """

    # Define constraint
    if sense=='>': 
        c = Constraint(x, lb=b)
    elif sense=='<': 
        c = Constraint(x, ub=b)
    else: 
        c = Constraint(x, lb=b, ub=b)
    
    # Return output
    return c

def cfr_optimize(cobra_model, on_list:list=[], off_list:list=[], 
                 on_params:set=(0.01, 0.001), off_params:set=(0.01, 0.), 
                 pfba_flag:bool=True, pfba_params:set=(1e-6, 0.), 
                 solver:str='gurobi', gurobi_params:dict=None, 
                 return_var:str='solution', verbose:bool=True): 
    """
    Optimize a COBRA model via constrain flux regulation (CFR). 
    
    This function simulates systems-level metabolism constrained by user-defined sets of active (on) and inactive (off) reactions. 
    
    Parameters
    ----------
    cobra_model : Model
        A COBRA Model object. 
    on_list : list
        List of active (on) genes or reactions.
    off_list : list
        List of inactive (off) genes or reactions.
    on_params : set, optional 
        Set of parameter constrains on active (on) reactions. The first argument corresponds to rho (weight coefficient) 
        while the second argument corresponds to epsilon 1 (minimum flux). 
    off_params : set, optional 
        Set of parameter constrains on inactive (off) reactions. The first argument corresponds to kappa (weight coefficient) 
        while the second argument corresponds to epsilon 2 (minimum flux). 
    pfba_flag : boolean, optional 
        Boolean flag for applying parsimonious flux balance analysis (pFBA).
    pfba_params : set, optional 
        Set of parameter constrains on non-inactive (non-off) reactions. The first argument corresponds to kappa2 (weight coefficient) 
        while the second argument corresponds to epsilon 3 (minimum flux). 
    solver : string, optional 
        String that specifies which optimization solver to use (only 'gurobi' and 'glpk' supported). 
    gurobi_params : dictionary, optional
        Dictionary of user-defined Gurobi model parameters.
    return_var : string, optional 
        String that specifies which variable to return (choose 'solution' or 'model'). 

    Returns
    -------
    Solution
        A COBRA solution object.
    Model
        A COBRA model object with CFR constraints applied.
        
    Raises
    ------
    TypeError
        Raised when non-float values are assigned to `on_params` or `off_params`.

    ValueError
        Raised when an incompatible value is provided as a function input.
        
    See Also
    --------
    `apply_cfr`: applies `cfr_optimize` across a set of conditions.
    
    Examples
    --------
    Usage cases of the `cfr_optimize` function.
    
    Case 1 (reaction lists): 
    >>> cobra_model = load_model('textbook')
    >>> on_list = cobra_model.reactions[:10]
    >>> off_list = cobra_model.reactions[20:30]
    >>> cfr_optimize(cobra_model, on_list, off_list)
    Optimal solution with objective value 0.874

    Case 2 (gene lists): 
    >>> cobra_model = load_model('textbook')
    >>> on_list = cobra_model.genes[:10]
    >>> off_list = cobra_model.genes[20:30]
    >>> cfr_optimize(cobra_model, on_list, off_list)
    Optimal solution with objective value 0.319

    Case 3 (custom constraint parameters): 
    >>> cobra_model = load_model('textbook')
    >>> on_list = cobra_model.genes[:10]
    >>> off_list = cobra_model.genes[20:30]
    >>> on_params = (0.1, 0.01)
    >>> off_params = (0.1, 0.01)
    >>> cfr_optimize(cobra_model, on_list, off_list, on_params=on_params, off_params=off_params)
    Optimal solution with objective value 0.197

    Case 4 (without pFBA): 
    >>> cobra_model = load_model('textbook')
    >>> on_list = cobra_model.genes[:10]
    >>> off_list = cobra_model.genes[20:30]
    >>> pfba_flag = False
    >>> cfr_optimize(cobra_model, on_list, off_list, pfba_flag=pfba_flag)
    Optimal solution with objective value 0.319

    Case 5 (using GLPK as solver): 
    >>> cobra_model = load_model('textbook')
    >>> on_list = cobra_model.genes[:10]
    >>> off_list = cobra_model.genes[20:30]
    >>> solver = 'glpk'
    >>> cfr_optimize(cobra_model, on_list, off_list, solver=solver)
    Optimal solution with objective value 0.324
    
    """
    # Check cobra_model
    cobra_attrs = ('genes', 'reactions', 'metabolites')
    if not all(hasattr(cobra_model, attr) for attr in cobra_attrs):
        raise ValueError("Invalid COBRA model: requires {} attributes.".format(', '.join(cobra_attrs)))

    # Check on_list
    if len(on_list)==0: 
        on_rxns = []
    else: 
        if any(item in cobra_model.reactions for item in on_list): 
            try: 
                on_rxns = [rxn.id for rxn in on_list]
            except: 
                on_rxns = on_list
        elif any(item in cobra_model.genes for item in on_list): 
            try:
                with cobra_model:
                    on_rxns = [rxn.id for rxn in geneKO_from_rxnGeneMat(cobra_model, on_list)]
                if verbose:
                    print('Determined active reactions from rxnGeneMat')
            except:
                with cobra_model: 
                    on_rxns = [rxn.id for rxn in knock_out_model_genes(cobra_model, on_list)]

    # Check off_list
    if len(off_list)==0: 
        off_rxns = []
    else: 
        if any(item in cobra_model.reactions for item in off_list): 
            try: 
                off_rxns = [rxn.id for rxn in off_list]
            except: 
                off_rxns = off_list
        elif any(item in cobra_model.genes for item in off_list): 
            try:
                with cobra_model:
                    off_rxns = [rxn.id for rxn in geneKO_from_rxnGeneMat(cobra_model, off_list)]
                if verbose:
                    print('Determined inactive reactions from rxnGeneMat')
            except:
                with cobra_model: 
                    off_rxns = [rxn.id for rxn in knock_out_model_genes(cobra_model, off_list)]

    # Consolidate conflicting on/off reactions (if needed)
    if len(set(on_rxns) & set(off_rxns)) > 0:
        off_rxns = [rxn for rxn in off_rxns if rxn not in on_rxns]
    on_rxns, off_rxns = sorted(on_rxns), sorted(off_rxns)

    # Check on_params
    if any(type(x) is not float for x in on_params): 
        raise TypeError('Provide a set of float values for on_params')
    if on_params[0] < 0 or on_params[0] > 10: 
        raise ValueError('Provide a float value between 0 and 10 for rho')
    if on_params[1] < 0 or on_params[1] > 1: 
        raise ValueError('Provide a float value between 0 and 1 for epsilon 1')
    
    # Check off_params
    if any(type(x) is not float for x in off_params): 
        raise TypeError('Provide a set of float values for off_params')
    if off_params[0] < 0 or off_params[0] > 10: 
        raise ValueError('Provide a float value between 0 and 10 for kappa')
    if off_params[1] < 0 or off_params[1] > 1: 
        raise ValueError('Provide a float value between 0 and 1 for epsilon 2')

    # Check pfba_params
    if any(type(x) is not float for x in pfba_params): 
        raise TypeError('Provide a set of float values for pfba_params')
    if pfba_params[0] < 0 or pfba_params[0] > 10: 
        raise ValueError('Provide a float value between 0 and 10 for kappa2')
    if pfba_params[1] < 0 or pfba_params[1] > 1: 
        raise ValueError('Provide a float value between 0 and 1 for epsilon 3')
    
    # Check solver
    if solver not in ('gurobi', 'glpk'): 
        raise ValueError('Invalid solver: must be either gurobi or glpk')

    # Check on gurobi_params
    if isinstance(gurobi_params, dict):
        raise ValueError('Provide a dictionary for gurobi_params')

    # Check return_var
    if return_var not in ('solution', 'model'): 
        raise ValueError('Invalid return_var: must be either solution or model')

    # Return default solution if both lists are empty
    if len(on_rxns)==0 and len(off_rxns)==0:
        if verbose:
            print('No CFR constraints detected: returning default solution')
        with cobra_model as model: 
            obj = model.solver.objective.expression
            s = str(obj).split(' ')
            obj_rxns = [i.split('*')[-1] for i in s if '*' in i]
            obj_rxns = [r for r in obj_rxns if r in model.reactions._dict.keys()]
            if pfba_flag: 
                variables = chain(*((rxn.forward_variable, rxn.reverse_variable) for rxn in model.reactions if rxn.id not in obj_rxns))
                model.objective.set_linear_coefficients({v: -w3 for v in variables})
            solution = model.optimize(solver)
            solution.objective_value = solution.fluxes[obj_rxns].sum()
    # Apply CFR
    else: 
        # Extract COBRA model data
        S = sp.csr_matrix(create_stoichiometric_matrix(cobra_model))
        r_dict = cobra_model.reactions._dict
        lb = np.array([rxn.lower_bound for rxn in cobra_model.reactions])
        ub = np.array([rxn.upper_bound for rxn in cobra_model.reactions])
        c = np.array([rxn.objective_coefficient for rxn in cobra_model.reactions])
        b = np.zeros(S.shape[0])

        # Define function inputs
        n1, n2 = len(on_rxns), len(off_rxns)
        w1, e1 = on_params
        w2, e2 = off_params
        dtype = float

        # Account for pFBA
        if pfba_flag: 
            pfba_rxns = sorted([rxn.id for rxn in cobra_model.reactions if rxn.id not in off_rxns])
        else: 
            pfba_rxns = []
        n3 = len(pfba_rxns)
        w3, e3 = pfba_params

        # Re-define inputs
        n = 2*n1 + 2*n2 + 2*n3
        r_list = np.repeat(on_rxns + off_rxns + pfba_rxns, 2)

        # Define CFR field: A matrix
        A0 = sp.csr_matrix((S.shape[0], n), dtype=dtype)
        row = np.arange(start=0, stop=n)
        col = np.array([r_dict[rxn] for rxn in r_list])
        val = np.ones(n, dtype=dtype)
        A1 = sp.csr_matrix((val, (row, col)), shape=(n, S.shape[1]))
        v1 = np.array([-(e1+10000), e1+10000], dtype=dtype)
        v2 = np.array([1., -1.], dtype=dtype)
        val = np.concatenate((np.tile(v1, n1), np.tile(v2, n2), np.tile(v2, n3)), dtype=dtype)
        A2 = sp.csr_matrix((val, (row, row)), shape=(n, n))
        A = sp.vstack([sp.hstack([S, A0]), sp.hstack([A1, A2])])

        # Define CFR field: rhs
        b1 = np.array([-10000, 10000], dtype=dtype)
        b2 = np.array([-e2, e2], dtype=dtype)
        b3 = np.array([-e3, e3], dtype=dtype)
        rhs = np.concatenate((b, np.tile(b1, n1), np.tile(b2, n2), np.tile(b3, n3)))

        # Define remaining CFR fields
        lb = np.concatenate((lb, np.zeros(n, dtype=dtype)))
        ub = np.concatenate((ub, np.ones(2*n1, dtype=dtype), 1000*np.ones(2*n2 + 2*n3, dtype=dtype)))
        vtype = np.repeat(['C', 'B', 'C'], (S.shape[1], 2*n1, 2*n2 + 2*n3))
        sense = np.concatenate((np.repeat(['='], S.shape[0]), np.tile(['>', '<'], n1 + n2 + n3)))
        obj = np.concatenate((c, np.repeat([w1, -w2, -w3], (2*n1, 2*n2, 2*n3))))

        # Construct CFR model
        if solver=='gurobi': 
            model = gp.Model(cobra_model.id + '_CFR')
            x = model.addMVar(shape=A.shape[1], lb=lb, ub=ub, vtype=vtype)
            model.addMConstr(A=A, x=x, sense=sense, b=rhs)
            model.setObjective(obj @ x, GRB.MAXIMIZE)
            if gurobi_params:
                for key, value in gurobi_params.items():
                    try:
                        model.setParam(key, value)
                    except GurobiError:
                        print(f'{key} is not a valid Gurobi parameter')
            model.Params.OutputFlag = 0
            model.optimize()
            if model.Status==2: 
                status = 'optimal'
                fluxes = pd.Series(model.X[:len(r_dict)], index=r_dict.keys())
                objective = fluxes.dot(c)
            else: 
                warn('Unable to determine optimal CFR solution. Returning indeterminate solution')
                status = 'not optimal'
                fluxes = pd.Series(np.nan, index=r_dict.keys())
                objective = np.nan
        elif solver=='glpk': 
            on_names = list(chain(*[('x_' + rxn + '_lower_on', 'x_' + rxn + '_upper_on') for rxn in on_rxns]))
            off_names = list(chain(*[('x_' + rxn + '_lower_off', 'x_' + rxn + '_upper_off') for rxn in off_rxns]))
            pfba_names = list(chain(*[('x_' + rxn + '_lower_pfba', 'x_' + rxn + '_upper_pfba') for rxn in pfba_rxns]))
            x_names = list(r_dict.keys()) + on_names + off_names + pfba_names
            x_list = [opt_variable(n, l, u, v) for n, l, u, v in tqdm(zip(x_names, lb, ub, vtype), desc='Defining variables', total=len(x_names), leave=False)]
            c_list = [opt_constraint(x_list[0], s, b) for s, b in tqdm(zip(sense, rhs), desc='Defining constraints', total=A.shape[0], leave=False)]
            model = Model('CFR')
            model.add(x_list + c_list)
            for i, constraint in tqdm(enumerate(model.constraints), desc='Updating constraints', total=len(model.constraints), leave=False): 
                idx = np.argwhere(A[i, :])[:, 1]
                if 0 not in idx: 
                    idx = np.append(idx, 0)
                x = [model.variables[i] for i in idx]
                c_dict = {v: w for v, w in zip(x, A[i, idx].toarray().reshape((-1,)))}
                constraint.set_linear_coefficients(c_dict)
            j = np.argwhere(c).reshape((-1,))[0]
            model.objective = Objective(c[j] * model.variables[j], direction='max')
            idx = np.arange(start=S.shape[1], stop=A.shape[1], step=1)
            x = [model.variables[i] for i in idx]
            o_dict = {v: w for v, w in zip(x, obj[idx])}
            model.objective.set_linear_coefficients(o_dict)
            model.configuration.timeout = 100
            status = model.optimize()
            if status=='optimal': 
                fluxes = pd.Series([model.variables[rxn].primal for rxn in r_dict.keys()], index=r_dict.keys())
                objective = fluxes.dot(c)
            else: 
                warn('Unable to determine optimal CFR solution. Returning indeterminate solution')
                fluxes = pd.Series(np.nan, index=r_dict.keys())
                objective = np.nan

        # Create COBRA solution
        solution = Solution(objective_value=objective, status=status, fluxes=fluxes)

    # Return output
    if return_var=='solution':
        return solution
    elif return_var=='model':
        return model

def apply_cfr(cobra_model, data:pd.DataFrame, thresh:set=(-2, 2), **kwargs): 
    """
    Apply CFR across a set of conditions.

    Parameters
    ----------
    cobra_model : Model
        A COBRA Model object.
    data : pandas.DataFrame
        A dataframe of omics data to infer active (on) and inactive (off) genes.
        None: the dataframe index is expected to match with genes defined in the COBRA model.
    thresh : set, optional
        A set of numeric values used to determine up- and down-regulated genes.
    kwargs : dict
        Keyword arguments accepted in the `cfr_optimize` function.

    Returns
    -------
    List
        A list containing Result objects for each condition.

    Raises
    ------
    ValueError
        Raised when an incompatible value is provided as a function input.
        
    See Also
    --------
    `cfr_optimize`: the function applied for each condition.
    `Result`: a class object containing input and output values for each condition.
    
    Examples
    --------
    Usage case of the `apply_cfr` function.

    >>> from random import choices
    >>> from cobra.io import load_model
    >>> cobra_model = load_model('textbook')
    >>> sample, k = range(-5, 5), len(cobra_model.genes)
    >>> data = pd.DataFrame({'C1': choices(sample, k=k), 'C2': choices(sample, k=k), 'C3': choices(sample, k=k)}, 
                            index=cobra_model.genes._dict.keys())
    >>> results = apply_cfr(cobra_model, data)
    >>> print(*['Objective for {} = {:.3f}\t'.format(key, value.solution.objective_value) for key, value in results.items()])
    Objective for C1 = 0.308	 Objective for C2 = 0.000	 Objective for C3 = 0.000

    """
    # Check cobra_model
    if not hasattr(cobra_model, 'genes'):
        raise ValueError('Invalid COBRA model: requires `genes` attribute')
    
    # Check thresh
    if any(not isinstance(x, (int, float)) for x in thresh): 
        raise TypeError('Provide a set of numeric values for `thresh`')
    
    # Extract model genes
    genes = [gene.id for gene in cobra_model.genes]

    # Filter data
    if not any(gene in genes for gene in data.index): 
        warn('Index of input data does not match COBRA model gene names')
    data = data[data.index.isin(genes)]

    # Determine up- and down-regulated genes
    off_dict = {col: data[data[col] < thresh[0]].index.tolist() for col in data.columns}
    on_dict = {col: data[data[col] > thresh[1]].index.tolist() for col in data.columns}

    # Determine solutions
    results = {}
    for col, up, down in tqdm(zip(data.columns, on_dict.values(), off_dict.values()), total=data.shape[1],
                              desc=f'Applying CFR across {data.shape[1]} conditions', leave=False):
        with cobra_model as model: 
            results[col] = Result(on_list=up, off_list=down, solution=cfr_optimize(model, up, down, **kwargs))

    # Return output
    return results

def summarize(results): 
    """
    Summarize output from `apply_cfr` into a dataframe.

    Parameters
    ----------
    results : list
        A list of Result objects (e.g., output from `apply_cfr`).

    Returns
    -------
    pandas.DataFrame
        A dataframe summarizing the number of up-regulated genes (N_up), 
        down-regulated genes (N_down), and objective value (Objective) for each condition in `results`.
    
    """
    # Create summary dataframe
    summary = pd.DataFrame({'N_up': [len(value.on_list) for value in results.values()], 
                            'N_down': [len(value.off_list) for value in results.values()], 
                            'Objective': [value.solution.objective_value for value in results.values()]}, 
                            index=results.keys())
    
    # Return output
    return summary

def get_fluxes(results, metadata:bool=False, cobra_model=None): 
    """
    Extract flux solutions for all conditions in the output for `apply_cfr`.

    Parameters
    ----------
    results : list
        A list of Result objects (e.g., output from `apply_cfr`).
    metadata : bool, optional
        Boolean flag to return metadata for reaction fluxes. 
    cobra_model : The COBRA Model object from which `results` was determined from.

    Returns
    -------
    pandas.DataFrame
        A dataframe of flux solutions for each condition in `results`.
    
    """
    # Extract flux solutions
    fluxes = {key: value.solution.fluxes for key, value in results.items()}

    # Create dataframe of fluxes
    df = pd.DataFrame(fluxes, index=next(iter(fluxes.values())).index)

    # Add metadata (if prompted)
    if metadata: 
        dfm = pd.DataFrame({'Reaction': [cobra_model.reactions.get_by_id(rxn).name for rxn in df.index]}, index=df.index)
        if all('subsystem' in vars(rxn).keys() for rxn in cobra_model.reactions): 
            dfm['Subsystem'] = [cobra_model.reactions.get_by_id(rxn).name for rxn in dfm.index]
        else: 
            dfm['Subsystem'] = np.nan
        df = dfm.merge(df, left_index=True, right_index=True)

    # Return output
    return df

def process_flux(cobra_model, fluxes, threshold=2): 
    """
    Process flux data into binarized format based on the reference COBRA model.

    Parameters
    ----------
    cobra_model : Model
        The COBRA Model object from which `fluxes` was determined from.
    fluxes : pandas.DataFrame
        A dataframe of flux solutions for a given set of conditions (e.g., output from `process_flux`).
    thresh : numeric (int or float)
        A numeric value used to identify positive entries in the binarized flux data.

    Returns
    -------
    pandas.DataFrame
        A dataframe of binarized flux solutions for each condition in `results`.
    
    """
    # Check cobra_model
    if not hasattr(cobra_model, 'reactions'): 
        raise ValueError('Invalid COBRA model: does not contain `reactions` attribute')
    
    # Check fluxes
    if not all(rxn1==rxn2 for rxn1, rxn2 in zip(cobra_model.reactions._dict.keys(), fluxes.index)): 
        raise ValueError('Index of `fluxes` does not match `cobra_model` reactions')
    
    # Check thresh
    if not isinstance(threshold, (int, float)): 
        raise ValueError('Provide a numeric value for `threshold`')

    # Add wild type flux solution
    df = fluxes.select_dtypes(include='number').copy()
    with HideOutput(): 
        df['WT'] = cfr_optimize(cobra_model).fluxes

    # Replace 0 with smallest positive flux
    df = df.replace(0, df[df > 0].min(axis=None))

    # Normalize fluxes based on WT
    dfz = df.iloc[:, :-1].sub(df['WT'], axis=0).div(df['WT'], axis=0)

    # Binarize based on threshold
    dfb = (dfz.abs() > threshold).astype(int)

    # Add metadata (if it exists)
    if fluxes.shape[1] > df.shape[1]: 
        dfm = fluxes.select_dtypes(exclude='number').copy()
        dfb = dfm.merge(dfb, left_index=True, right_index=True)

    # Return output
    return dfb

class Result:
    """
    A unified interface to access the results from `apply_cfr`.

    Parameters
    ----------
    on_list : list
        Contains the list of up-regulated genes.
    off_list : list
        Contains the list of down-regulated genes.
    solution : COBRA Solution
        The COBRA solution after applying CFR.

    Attributes
    ----------
    on_list : list
        Contains the list of up-regulated genes.
    off_list : list
        Contains the list of down-regulated genes.
    solution : COBRA Solution
        The COBRA solution after applying CFR.

    Notes
    -----
    Result is meant to be constructed by `apply_cfr`. Please look at that
    function to fully understand the `Result` class.

    """
    def __init__(
        self,
        on_list: list,
        off_list: list,
        solution : Solution,
    ) -> None:
        """
        Initialize a `Result` from its components.

        """
        super().__init__()
        self.on_list = on_list
        self.off_list = off_list
        self.solution = solution

class HideOutput:
    """
    Hide print statements for a given function call.
    
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout