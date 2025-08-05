"""
This module provides methods for the Python implementation of the constrain flux regulation (CFR) algorithm.

"""

__licence__ = "GPL GNU"
__docformat__ = "reStructuredText"

# Dependencies
import os, sys
import optlang
import numpy as np
import pandas as pd
from tqdm import tqdm
from sympy.matrices import SparseMatrix
import gurobipy as gp
from gurobipy import GRB
import scipy.sparse as sp
from itertools import chain
from warnings import warn
from importlib import import_module
from cobra.core import Solution, get_solution
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.util.array import create_stoichiometric_matrix
from cobra.manipulation.delete import knock_out_model_genes
from optlang.glpk_interface import Variable, Constraint, Objective, Model

# Methods
def cfr(cobra_model, on_rxns:list=None, off_rxns:list=None, 
        on_gns:list=None, off_gns:list=None, 
        on_params:set=(0.01, 0.001), off_params:set=(0.01, 0.001), 
        pfba_flag:bool=True, solver:str=None, silent:bool=False): 
    
    """Function for constrain flux regulation (CFR). 
    
    This function simulates systems-level metabolism constrained by user-defined sets of active (on) and inactive (off) reactions. 
    
    Parameters
    ----------
    cobra_model : DictList
        Description of arg1. 
    on_rxns : list, optional
        List of active (on) reactions.
    off_rxns : list
        List of inactive (off) reactions.
    on_gns : list
        List of up-regulated genes, associated with active (on) reactions. 
    off_gns : list
        List of down-regulated genes, associated with inactive (off) reactions.
    on_params : set
        Set of parameter constrains on active reactions. The first argument corresponds to rho (weight coefficient) while the second argument corresponds to epsilon 1 (minimum flux). 
    off_params : set
        Set of parameter constrains on inactive reactions. The first argument corresponds to kappa (weight coefficient) while the second argument corresponds to epsilon 2 (minimum flux). 
    pfba_flag : boolean
        Boolean flag for applying parsimonious flux balance analysis (pFBA).
    solver : string, optional
        Boolean flag for applying flux variability analysis (FVA). 
        
    Returns
    -------
    type
        COBRA solution object.
        
    Raises
    ------
    ErrorType
        Error description. 
        
    See Also
    --------
    otherfunc: other related function.
    
    Examples
    --------
    Usage cases of the `cfr` function.
    
    >>> line 1
    >>> line 2
    Return value
    
    """
    # Check inputs
    attrs = vars(cobra_model)
    if not any(var in attrs.keys() for var in ('genes', 'reactions', 'metabolites')): 
        raise AssertionError('Provide a valid COBRA model with gene, reaction, and metabolite attributes')
    elif all(var in attrs.keys() for var in ('genes', 'reactions', 'metabolites')): 
        model = cobra_model.copy()
        if silent is False: 
            print('Valid COBRA model provided')

    if on_rxns is None: 
        if on_gns is None: 
            raise ValueError('Provide a list of active reactions or up-regulated genes')
        elif not any(gene in model.genes for gene in on_gns): 
            warn('None of the up-regulated genes present in the given model')
            on_rxns = []
        else: 
            if not all(gene in model.genes for gene in on_gns): 
                warn('Not all up-regulated genes are present in the given model')
            elif all(gene in model.genes for gene in on_gns): 
                if silent is False: 
                    print('All up-regulated genes accounted in given model')
            with model: 
                on_rxns = knock_out_model_genes(model, on_gns)
            # on_rxns = list(set().union(*[model.genes.get_by_id(g).reactions for g in on_gns]))
            if silent is False: 
                print('Active reactions inferred from up-regulated genes')
    elif not any(rxn in model.reactions for rxn in on_rxns): 
        raise AssertionError('None of the active reactions present in the given model')
    elif not all(rxn in model.reactions for rxn in on_rxns): 
        warn('Not all reactions in active list are present in the given model')
    elif all(rxn in model.reactions for rxn in on_rxns): 
        if silent is False: 
            print('All active reactions accounted in given model')
    
    if off_rxns is None: 
        if off_gns is None: 
            raise ValueError('Provide a list of inactive reactions or down-regulated genes')
        elif not any(gene in model.genes for gene in off_gns): 
            warn('None of the down-regulated genes present in the given model')
            off_rxns = []
        else: 
            if not all(gene in model.genes for gene in off_gns): 
                warn('Not all down-regulated genes are present in the given model')
            elif all(gene in model.genes for gene in off_gns): 
                if silent is False: 
                    print('All down-regulated genes accounted in given model')
            with model: 
                off_rxns = knock_out_model_genes(model, off_gns)
            # off_rxns = list(set().union(*[model.genes.get_by_id(g).reactions for g in off_gns]))
            if silent is False: 
                print('Inactive reactions inferred from down-regulated genes')
    elif not any(rxn in model.reactions for rxn in off_rxns): 
        raise AssertionError('None of the inactive reactions present in the given model')
    elif not all(rxn in model.reactions for rxn in off_rxns): 
        warn('Not all reactions in inactive list are present in the given model')
    elif all(rxn in model.reactions for rxn in off_rxns): 
        if silent is False: 
            print('All inactive reactions accounted in given model')

    if any(type(x) is not float for x in on_params): 
        raise TypeError('Provide a set of float values for on_params')
    if on_params[0] < 0 or on_params[0] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for rho')
    if on_params[1] < 0 or on_params[1] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for epsilon 1')
    
    if any(type(x) is not float for x in off_params): 
        raise TypeError('Provide a set of float values for off_params')
    if off_params[0] < 0 or off_params[0] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for kappa')
    if off_params[1] < 0 or off_params[1] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for epsilon 2')

    if solver is None: 
        from optlang import Variable, Constraint, Objective
    else: 
        solver_dict = optlang.available_solvers
        solvers = [key.lower() for key, value in solver_dict.items() if value is True]
        if solver not in solvers: 
            raise ValueError('Provide a valid and available solver\nAvailable solvers: {}'.format(', '.join(solvers)))
        else: 
            model.solver = solver
            solver_mod = import_module('optlang.{}_interface'.format(solver))
            Variable = getattr(solver_mod, 'Variable')
            Constraint = getattr(solver_mod, 'Constraint')
            Objective = getattr(solver_mod, 'Objective')
            if silent is False: 
                print('{} methods successfully loaded'.format(solver.upper()))
    
    # Extract objective information
    obj = model.solver.objective.expression
    s = str(obj).split(' ')
    obj_rxns = [i.split('*')[-1] for i in s if '*' in i]
    obj_rxns = [r for r in obj_rxns if r in model.reactions._dict.keys()]
    print(obj_rxns)

    # Instantiate function variables
    i, M, x, cons = 0, 10000, list(), list()

    # Add constraints for active (on) reactions
    r, e = on_params
    for rxn in tqdm(on_rxns, desc='Adding constraints for active (on) reactions'): 
        # Get reaction variable
        v = model.solver.variables[rxn.id]
        # Create lower constraint
        x.append(Variable('x' + str(i), lb=0, ub=1, type='binary'))
        cons.append(Constraint(1*v - (e + M)*x[i], lb=-M))
        obj = obj + r*x[i]
        i = i + 1
        # Create upper constraint
        x.append(Variable('x' + str(i), lb=0, ub=1, type='binary'))
        cons.append(Constraint(1*v + (e + M)*x[i], ub=M))
        obj = obj + r*x[i]
        i = i + 1

    # Add constraints for inactive (off) reactions
    k, e = off_params
    for rxn in tqdm(off_rxns, desc='Adding constraints for inactive (off) reactions'): 
        # Get reaction variable
        v = model.solver.variables[rxn.id]
        # Create lower constraint
        x.append(Variable('x' + str(i), lb=0, ub=1000, type='continuous'))
        cons.append(Constraint(1*v + 1*x[i], lb=-e))
        obj = obj - k*x[i]
        i = i + 1
        # Create upper constraint
        x.append(Variable('x' + str(i), lb=0, ub=1000, type='continuous'))
        cons.append(Constraint(1*v - 1*x[i], ub=e))
        obj = obj - k*x[i]
        i = i + 1

    """# Add pFBA constraints (if prompted)
    if pfba_flag is True: 
        if silent is False: 
            print('Applying pFBA constraints')
        pfba_rxns = [rxn for rxn in model.reactions if rxn not in off_rxns]
        k, e = 1e-6, 0
        for rxn in tqdm(pfba_rxns, desc='Adding constraints for pFBA'): 
            # Get reaction variable
            v = model.solver.variables[rxn._id]
            # Create lower constraint
            x.append(Variable('x' + str(i), lb=0, ub=1000, type='continuous'))
            cons.append(Constraint(1*v + 1*x[i], lb=-e))
            obj = obj - k*x[i]
            i = i + 1
            # Create upper constraint
            x.append(Variable('x' + str(i), lb=0, ub=1000, type='continuous'))
            cons.append(Constraint(1*v - 1*x[i], ub=e))
            obj = obj - k*x[i]
            i = i + 1"""

    # Update model
    model.solver.add(cons)
    # model.solver.objective = Objective(obj, direction='max')

    # Add pFBA constraint (if prompted)
    if pfba_flag is True: 
        add_pfba(model, objective=Objective(obj, direction='max'))
    else: 
        model.solver.objective = Objective(obj, direction='max')
    print(obj)

    # Print updated model specs (if prompted)
    if silent is False: 
        print('Updated model specs:')
        print('\tNo. of genes: {}'.format(len(model.genes)))
        print('\tNo. of reactions: {}'.format(len(model.reactions)))
        print('\tNo. of metabolites: {}'.format(len(model.metabolites)))
        print('\tNo. of solver variables: {}'.format(len(model.solver.variables)))
        print('\tNo. of solver constraints: {}'.format(len(model.solver.constraints)))

    # Determine solution
    model.optimize()

    # Return solution according to original model
    solution = get_solution(model)
    solution.objective_value = solution.fluxes[obj_rxns].sum()
    return solution

def cfr_variables(model, rxn, status): 
    """
    Defines CFR variables for a reaction ID and its status.
    """
    # Define variable name
    try: 
        v_name = '_'.join(['x', rxn.id, status])
    except: 
        v_name = '_'.join(['x', rxn, status])

    # Define variables
    if status=='on': 
        x_lower = model.problem.Variable(v_name + '_lower', lb=0, ub=1, type='binary')
        x_upper = model.problem.Variable(v_name + '_upper', lb=0, ub=1, type='binary')
    else: 
        x_lower = model.problem.Variable(v_name + '_lower', lb=0, ub=1000, type='continuous')
        x_upper = model.problem.Variable(v_name + '_upper', lb=0, ub=1000, type='continuous')

    # Return output as tuple
    return (x_lower, x_upper)

def cfr_constraints(model, rxn, cfr_vars, epsilon): 
    """
    Defines CFR constraints for a given CFR variable.
    
    """
    # Find rxn variable
    try: 
        v = model.variables[rxn.id]
    except: 
        v = model.variables[rxn]

    # Resolve CFR variables
    x_lower, x_upper = cfr_vars

    # Define contraint
    if 'on' in x_lower.name: 
        c_lower = model.problem.Constraint(1*v - (epsilon + 10000)*x_lower, lb=-10000)
        c_upper = model.problem.Constraint(1*v + (epsilon + 10000)*x_upper, ub=10000)
    else: 
        c_lower = model.problem.Constraint(1*v + 1*x_lower, lb=-epsilon)
        c_upper = model.problem.Constraint(1*v - 1*x_upper, ub=epsilon)
    
    # Return output as tuple
    return (c_lower, c_upper)

def cfr_optimize_deprecated(model, on_list:list=None, off_list:list=None, 
        on_params:set=(0.01, 0.001), off_params:set=(0.01, 0.001), 
        pfba_flag:bool=True, solver:str=None): 
    
    """
    Function for constrain flux regulation (CFR). 
    
    This function simulates systems-level metabolism constrained by user-defined sets of active (on) and inactive (off) reactions. 
    
    Parameters
    ----------
    model : DictList
        Description of arg1. 
    on_list : list
        List of active (on) genes or reactions.
    off_list : list
        List of inactive (off) genes or reactions.
    on_params : set, optional (Default: on_params=(0.01, 0.001))
        Set of parameter constrains on active reactions. The first argument corresponds to rho (weight coefficient) while the second argument corresponds to epsilon 1 (minimum flux). 
    off_params : set, optional (Default: off_params=(0.01, 0.001))
        Set of parameter constrains on inactive reactions. The first argument corresponds to kappa (weight coefficient) while the second argument corresponds to epsilon 2 (minimum flux). 
    pfba_flag : boolean, optional (Default: pfba_flag=True)
        Boolean flag for applying parsimonious flux balance analysis (pFBA).
    solver : string, optional (Default: solver=None)
        String that specifies which optimization solver to use. 
        
    Returns
    -------
    Solution
        COBRA solution object.
        
    Raises
    ------
    ErrorType
        Error description. 
        
    See Also
    --------
    otherfunc: other related function.
    
    Examples
    --------
    Usage cases of the `cfr` function.
    
    >>> line 1
    >>> line 2
    Return values
    
    """
    # Check cobra_model
    cobra_attrs = ('problem', 'solver', 'genes', 'reactions', 'metabolites')
    if not all(hasattr(model, attr) for attr in cobra_attrs):
        raise ValueError("Invalid COBRA model: Missing {} attribute.".format(', '.join(cobra_attrs)))
    '''with HideOutput(): 
        # model = cobra_model.copy()
        model = cobra_model'''
    if solver=='gurobi' or 'gurobi' in str(type(model.solver)): 
        model.solver.problem.Params.IterationLimit = 1e8
        model.solver.problem.Params.FeasibilityTol = 1e-6
        model.solver.problem.Params.IntFeasTol = 1e-5
        model.solver.problem.Params.OptimalityTol = 1e-6
        model.solver.problem.Params.Presolve = -1
        model.solver.problem.Params.Method = -1

    # Check on_list
    if on_list is None: 
        raise ValueError('Provide a list of active reactions or up-regulated genes')
    elif len(on_list)==0: 
        on_rxns = []
    else: 
        if any(item in model.reactions for item in on_list): 
            on_rxns = list(on_list)
        elif any(item in model.genes for item in on_list): 
            with model: 
                on_rxns = knock_out_model_genes(model, on_list)

    # Check off_list
    if off_list is None: 
        raise ValueError('Provide a list of inactive reactions or down-regulated genes')
    elif len(off_list)==0: 
        off_rxns = []
    else: 
        if any(item in model.reactions for item in off_list): 
            off_rxns = list(off_list)
        elif any(item in model.genes for item in off_list): 
            with model: 
                off_rxns = knock_out_model_genes(model, off_list)

    # Check on_params
    if any(type(x) is not float for x in on_params): 
        raise TypeError('Provide a set of float values for on_params')
    if on_params[0] < 0 or on_params[0] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for rho')
    if on_params[1] < 0 or on_params[1] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for epsilon 1')
    
    # Check off_params
    if any(type(x) is not float for x in off_params): 
        raise TypeError('Provide a set of float values for off_params')
    if off_params[0] < 0 or off_params[0] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for kappa')
    if off_params[1] < 0 or off_params[1] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for epsilon 2')
    
    # Check solver
    if solver is None: 
        from optlang import Variable, Constraint
    else: 
        solver_dict = optlang.available_solvers
        solvers = [key.lower() for key, value in solver_dict.items() if value is True]
        if solver not in solvers: 
            raise ValueError('Provide a valid and available solver\nAvailable solvers: {}'.format(', '.join(solvers)))
        else: 
            model.solver = solver
            solver_mod = import_module('optlang.{}_interface'.format(solver))
            Variable = getattr(solver_mod, 'Variable')
            Constraint = getattr(solver_mod, 'Constraint')

    # Return default solution if both lists are empty
    if len(on_rxns)==0 and len(off_rxns)==0: 
        warn('No CFR constraints detected. Returning default solution')
        if solver==None: 
            solution = model.optimize()
        else: 
            solution = model.optimize(solver)
    else: 
        # Extract objective information
        obj = model.objective.expression
        s = str(obj).split(' ')
        obj_rxns = [i.split('*')[-1] for i in s if '*' in i]
        obj_rxns = [r for r in obj_rxns if r in model.reactions._dict.keys()]

        # Define function inputs
        n1, n2 = len(on_rxns), len(off_rxns)
        w1, e1 = on_params
        w2, e2 = off_params
        status = ['on']*n1 + ['off']*n2
        weights = [w1]*2*n1 + [-w2]*2*n2
        epsilon = [e1]*n1 + [e2]*n2
        rxn_list = on_rxns + off_rxns

        # Account for pFBA 
        if pfba_flag: 
            pfba_rxns = [rxn for rxn in model.reactions if rxn not in off_rxns]
            status = status + ['pfba']*len(pfba_rxns)
            weights = weights + [-1e-6]*2*len(pfba_rxns)
            epsilon = epsilon + [0]*len(pfba_rxns)
            rxn_list = rxn_list + pfba_rxns

        # Define CFR variables
        x_list = [cfr_variables(model, rxn, state) for rxn, state in zip(rxn_list, status)]
        x_all = [x for x in chain(*x_list)]

        # Define CFR constraints
        c_list = [cfr_constraints(model, rxn, x, e) for rxn, x, e in zip(rxn_list, x_list, epsilon)]
        cons = [c for c in chain(*c_list)]

        # Add CFR constraints
        model.add_cons_vars(x_all + cons)

        # Update objective
        model.objective.set_linear_coefficients({v: w for v, w in zip(x_all, weights)})

        # Determine solution
        solution = model.optimize()

        # Return solution according to original model
        solution.objective_value = solution.fluxes[obj_rxns].sum()

    # Return output
    return solution

def opt_variable(name, lb, ub, vtype): 
    """
    Defines an optlang variable.
    """
    if vtype=='B': 
        x = Variable(name, lb=lb, ub=ub, type='binary')
    else: 
        x = Variable(name, lb=lb, ub=ub, type='continuous')

    return x

def opt_constraint(x, sense, b): 
    """
    Defines an optlang constraint.
    
    """
    # idx = np.argwhere(a)[:, 1]
    # expr = np.sum(a[0, idx].toarray() * x[idx])
    '''if len(idx) < 10: 
        expr = a[0, idx].toarray().dot(x[idx])[0]
    else: 
        expr = a[0, idx[0]].toarray().dot(x[idx[0]])[0]
    j = np.argwhere(a)[0, 1]
    expr = a[0, j] * x[j]'''
    if sense=='>': 
        c = Constraint(x, lb=b)
    elif sense=='<': 
        c = Constraint(x, ub=b)
    else: 
        c = Constraint(x, lb=b, ub=b)
    
    # Return output
    return c

def set_coefficients(c, a, x): 
    """
    Set linear coefficients to given constraint.
    
    """
    idx = np.argwhere(a)[:, 1]
    c_dict = {v: w for v, w in zip(x[idx], a[0, idx].toarray().reshape((-1,)))}
    c.set_linear_coefficients(c_dict)

    # Return output
    return c

def cfr_optimize(cobra_model, on_list:list=None, off_list:list=None, on_params:set=(0.01, 0.001), off_params:set=(0.01, 0.001), 
        pfba_flag:bool=True, solver:str='gurobi'): 
    
    """
    Function for constrain flux regulation (CFR). 
    
    This function simulates systems-level metabolism constrained by user-defined sets of active (on) and inactive (off) reactions. 
    
    Parameters
    ----------
    model : DictList
        Description of arg1. 
    on_list : list
        List of active (on) genes or reactions.
    off_list : list
        List of inactive (off) genes or reactions.
    on_params : set, optional (Default: on_params=(0.01, 0.001))
        Set of parameter constrains on active reactions. The first argument corresponds to rho (weight coefficient) while the second argument corresponds to epsilon 1 (minimum flux). 
    off_params : set, optional (Default: off_params=(0.01, 0.001))
        Set of parameter constrains on inactive reactions. The first argument corresponds to kappa (weight coefficient) while the second argument corresponds to epsilon 2 (minimum flux). 
    pfba_flag : boolean, optional (Default: pfba_flag=True)
        Boolean flag for applying parsimonious flux balance analysis (pFBA).
    solver : string, optional (Default: solver=None)
        String that specifies which optimization solver to use. 
        
    Returns
    -------
    Solution
        COBRA solution object.
        
    Raises
    ------
    ErrorType
        Error description. 
        
    See Also
    --------
    otherfunc: other related function.
    
    Examples
    --------
    Usage cases of the `cfr` function.
    
    >>> line 1
    >>> line 2
    Return values
    
    """
    # Check cobra_model
    cobra_attrs = ('genes', 'reactions', 'metabolites')
    if not all(hasattr(cobra_model, attr) for attr in cobra_attrs):
        raise ValueError("Invalid COBRA model: Missing {} attribute.".format(', '.join(cobra_attrs)))

    # Check on_list
    if on_list is None: 
        raise ValueError('Provide a list of active reactions or up-regulated genes')
    elif len(on_list)==0: 
        on_rxns = []
    else: 
        if any(item in cobra_model.reactions for item in on_list): 
            try: 
                on_rxns = [rxn.id for rxn in on_list]
            except: 
                on_rxns = on_list
        elif any(item in cobra_model.genes for item in on_list): 
            with cobra_model: 
                on_rxns = [rxn.id for rxn in knock_out_model_genes(cobra_model, on_list)]

    # Check off_list
    if off_list is None: 
        raise ValueError('Provide a list of inactive reactions or down-regulated genes')
    elif len(off_list)==0: 
        off_rxns = []
    else: 
        if any(item in cobra_model.reactions for item in off_list): 
            try: 
                off_rxns = [rxn.id for rxn in off_list]
            except: 
                off_rxns = off_list
        elif any(item in cobra_model.genes for item in off_list): 
            with cobra_model: 
                off_rxns = [rxn.id for rxn in knock_out_model_genes(cobra_model, off_list)]

    # Check on_params
    if any(type(x) is not float for x in on_params): 
        raise TypeError('Provide a set of float values for on_params')
    if on_params[0] < 0 or on_params[0] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for rho')
    if on_params[1] < 0 or on_params[1] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for epsilon 1')
    
    # Check off_params
    if any(type(x) is not float for x in off_params): 
        raise TypeError('Provide a set of float values for off_params')
    if off_params[0] < 0 or off_params[0] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for kappa')
    if off_params[1] < 0 or off_params[1] > 1: 
        raise ValueError('Provide a boolean value between 0 and 1 for epsilon 2')
    
    # Check solver
    if solver not in ('gurobi', 'glpk'): 
        raise ValueError('Invalid solver. Must be either gurobi or glpk')

    # Return default solution if both lists are empty
    if len(on_rxns)==0 and len(off_rxns)==0: 
        warn('No CFR constraints detected. Returning default solution')
        solution = cobra_model.optimize(solver)
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
            pfba_rxns = [rxn.id for rxn in cobra_model.reactions if rxn.id not in off_rxns]
        else: 
            pfba_rxns = []
        n3, w3, e3 = len(pfba_rxns), 1e-6, 0.

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
            model.Params.OutputFlag = 0
            model.optimize()
            if model.Status==2: 
                status = 'optimal'
                fluxes = pd.Series(model.X[:len(r_dict)], index=r_dict.keys())
                objective = fluxes[c!=0].sum()
            else: 
                warn('Unable to determine optimal CFR solution. Returning indeterminate solution')
                status = 'not optimal'
                fluxes = pd.Series(np.nan, index=r_dict.keys())
                objective = np.nan
        elif solver=='glpk': 
            on_names = list(chain(*[('x_' + rxn + '_lower', 'x_' + rxn + '_upper') for rxn in on_rxns]))
            off_names = list(chain(*[('x_' + rxn + '_lower', 'x_' + rxn + '_upper') for rxn in off_rxns]))
            pfba_names = list(chain(*[('x_' + rxn + '_lower_pfba', 'x_' + rxn + '_upper_pfba') for rxn in pfba_rxns]))
            x_names = list(r_dict.keys()) + on_names + off_names + pfba_names
            x_list = [opt_variable(n, l, u, v) for n, l, u, v in tqdm(zip(x_names, lb, ub, vtype), desc='Defining variables', total=len(x_names))]
            c_list = [opt_constraint(x_list[0], s, b) for s, b in tqdm(zip(sense, rhs), desc='Defining constraints', total=A.shape[0])]
            model = Model('CFR')
            model.add(x_list + c_list)
            for i, constraint in tqdm(enumerate(model.constraints), desc='Updating constraints', total=len(model.constraints)): 
                idx = np.argwhere(A[i, :])[:, 1]
                if 0 not in idx: 
                    idx = np.append(idx, 0)
                x = [model.variables[i] for i in idx]
                c_dict = {v: w for v, w in zip(x, A[i, idx].toarray().reshape((-1,)))}
                constraint.set_linear_coefficients(c_dict)
            print('Defining objective')
            j = np.argwhere(c).reshape((-1,))[0]
            model.objective = Objective(c[j] * model.variables[j], direction='max')
            idx = np.arange(start=S.shape[1], stop=A.shape[1], step=1)
            x = [model.variables[i] for i in idx]
            o_dict = {v: w for v, w in zip(x, obj[idx])}
            model.objective.set_linear_coefficients(o_dict)
            print('Objective defined')
            status = model.optimize()
            if status=='optimal': 
                fluxes = pd.Series([model.variables[rxn].primal for rxn in r_dict.keys()], index=r_dict.keys())
                objective = fluxes[c!=0].sum()
            else: 
                warn('Unable to determine optimal CFR solution. Returning indeterminate solution')
                fluxes = pd.Series(np.nan, index=r_dict.keys())
                objective = np.nan

        # Create COBRA solution
        solution = Solution(objective_value=objective, status=status, fluxes=fluxes)

    # Return output
    return solution

def apply_cfr(cobra_model, data:pd.DataFrame, thresh=(-2, 2), **kwargs): 
    """
    Apply CFR given a COBRA model and an omics dataframe.
    
    """
    # Extract model genes
    genes = [gene.id for gene in cobra_model.genes]

    # Filter data
    data = data[data.index.isin(genes)]

    # Determine up- and down-regulated genes
    off_dict = {col: data[data[col] < thresh[0]].index.tolist() for col in data.columns}
    on_dict = {col: data[data[col] > thresh[1]].index.tolist() for col in data.columns}

    # Determine solutions
    results = {}
    for col, up, down in tqdm(zip(data.columns, on_dict.values(), off_dict.values()), total=data.shape[1],
                              desc='Applying CFR across {} conditions'.format(data.shape[1])):
        with cobra_model as model: 
            results[col] = Result(on_list=up, off_list=down, solution=cfr_optimize(model, up, down, **kwargs))

    # Define output
    return results

def summarize(results): 
    """
    Summarize output from `apply_cfr` into a dataframe.
    
    """
    # Create summary dataframe
    summary = pd.DataFrame({'N_up': [len(value.on_list) for value in results.values()], 
                            'N_down': [len(value.off_list) for value in results.values()], 
                            'Objective': [value.solution.objective_value for value in results.values()]}, 
                            index=results.keys())
    
    # Return output
    return summary

def get_fluxes(results): 
    """
    Extract flux solutions for all conditions in the output for `apply_cfr`.
    
    """
    # Extract flux solutions
    fluxes = {key: value.solution.fluxes for key, value in results.items()}

    # Create dataframe of fluxes
    df = pd.DataFrame(fluxes, index=next(iter(fluxes.values())).index)

    # Return output
    return df

def process_flux(fluxes, cobra_model, threshold:float=2): 
    """
    Process flux data into binarized data based on a central COBRA model.
    
    """
    # Add wild type flux solution
    df = fluxes.copy()
    # df['WT'] = cobra_model.optimize().fluxes
    with HideOutput(): 
        df['WT'] = cfr_optimize(cobra_model, [], []).fluxes

    # Replace 0 with smallest positive flux
    df = df.replace(0, df[df > 0].min(axis=None))

    # Normalize fluxes based on WT
    dfz = df.iloc[:, :-1].sub(df['WT'], axis=0).div(df['WT'], axis=0)

    # Binarize based on threshold
    dfb = (dfz.abs() > threshold).astype(int)

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