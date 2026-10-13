function [solution, gurobiModel] = cfr_optimize(model, onList, offList, options)
    % CFR_OPTIMIZE Optimize a COBRA model via constrain flux regulation (CFR).
    %
    % DESCRIPTION
    % This function simulates flux through reactions annotated in a
    % genome-scale metabolic model (GEM) by constraining the model using
    % information on "on" and "off" reactions (can be based directly or by
    % differentially expressed gene (DEG) data.
    %
    % Default values match cfr_optimize in cfr_suite.py, and the outputs
    % cover both of its return_var choices (a solution and the solver
    % model). Differences from
    % constrain_flux_regulation.m are listed in cfr_optimize_legacy.m
    % (both in matlab/legacy); in
    % addition, this version uses name-value options, returns the final
    % reaction lists and per-reaction diagnostics, and adds a pFBA tolerance.
    %
    % STEPS
    % 1. Resolve gene/reaction lists and set parameters
    % 2. Convert model to gurobi format
    % 3. Add constraints for "on", "off" and pFBA reactions
    % 4. Define function outputs
    %
    % Author:   Sriram Chandrasekaran
    % Created:  October 12, 2010
    % Updates:  October 8, 2026 (Carolina H. Chung)
    %           October 6, 2026 (Carolina H. Chung)
    %           July 14, 2026 (Carolina H. Chung)
    %           September 20, 2020 (Carolina H. Chung)
    %           November 6, 2019 (Scott Campit)
    %           November 8, 2018 (Fangzhou Shen)

    % I/O
    %{
    REQUIRED INPUTS:
        1. model:       genome-scale metabolic model (GEM) with fields
                        S, rxns, lb, ub and c (b, csense and genes/rules
                        are used when present)
           -> type:     MATLAB structure

        2. onList:      "on" reactions or up-regulated genes
           -> type:     cell array or string array

        3. offList:     "off" reactions or down-regulated genes
           -> type:     cell array or string array

    OPTIONS (name-value pairs; Greek symbols as in the CFR literature):
        listType:       "auto" (detect from IDs, as in cfr_suite.py),
                        "reactions" or "genes"
           -> default:  "auto"
        onWeight:       weight for "on" reactions (rho), 0-10
           -> default:  0.01
        onMinFlux:      minimum |flux| for "on" reactions (epsilon 1), 0-1
           -> default:  1e-3
        offWeight:      weight for "off" reactions (kappa), 0-10
           -> default:  0.01
        offTolerance:   flux tolerance for "off" reactions (epsilon 2), 0-1
           -> default:  0
        pfba:           apply parsimonious FBA to all non-"off" reactions
           -> default:  true
        pfbaWeight:     weight for pFBA reactions (kappa 2), 0-10
           -> default:  1e-6
        pfbaTolerance:  flux tolerance for pFBA reactions (epsilon 3), 0-1
           -> default:  0
        gurobiParams:   gurobi parameters, validated with
                        cfr_suite.check_gurobi_params (invalid or duplicated
                        parameters raise a warning and are ignored);
                        OutputFlag is always set to 0; IntFeasTol
                        defaults to 1e-9 unless given
           -> default:  struct()
        verbose:        print progress statements and valid gurobi
                        parameters (as verbose in cfr_suite.py)
           -> default:  true
        checkGurobiParams:  validate gurobiParams; apply_cfr sets false
                        after validating them once for all conditions
           -> default:  true

        onWeight, onMinFlux, offWeight and offTolerance are scalars or
        vectors with one value per entry of a reaction list. With gene
        lists, use scalars.

    OUTPUTS:
        1. solution:        COBRA solution structure (fields as in
                            optimizeCbModel, plus CFR diagnostics);
                            return_var='solution' in cfr_suite.py
           -> f:            value of the model objective (model.c' * v);
                            objective_value in cfr_suite.py
           -> v:            flux solution (order of model.rxns); NaN if not
                            optimal; fluxes in cfr_suite.py
           -> x:            same as v (legacy COBRA field)
           -> stat:         1 optimal, 0 infeasible, 2 unbounded,
                            3 suboptimal, -1 other (time limit, etc.)
           -> origStat:     gurobi status (e.g. 'OPTIMAL', 'INFEASIBLE')
           -> solver:       'gurobi'
           -> cfrObjective: value of the full CFR objective
           -> onReactions:  final "on" reactions (after gene mapping,
                            removal of unknown IDs and duplicates)
           -> offReactions: final "off" reactions (without pFBA reactions)
           -> onSatisfied:  logical, |flux| >= onMinFlux per "on" reaction
           -> offSatisfied: logical, |flux| <= offTolerance per "off" reaction

        2. gurobiModel:     gurobi model structure with CFR constraints;
                            return_var='model' in cfr_suite.py. Solve it
                            again with gurobi(gurobiModel, struct(...
                            'OutputFlag', 0, 'IntFeasTol', 1e-9)).

    EXAMPLE USAGE:
        1. Default usage:
           >> solution = cfr_suite.cfr_optimize(model, onRxns, offRxns);
           >> solution.f, solution.v

        2. Gene lists with custom weights:
           >> solution = cfr_suite.cfr_optimize(model, upGenes, downGenes, ...
                    listType="genes", onWeight=1e-3, offWeight=1e-3)

        3. Without pFBA (name-value syntax for MATLAB R2019b-R2020b):
           >> solution = cfr_suite.cfr_optimize(model, onRxns, offRxns, 'pfba', false)

        4. Solution and gurobi model with CFR constraints:
           >> [solution, gurobiModel] = cfr_suite.cfr_optimize(model, onRxns, offRxns);
           >> [~, gurobiModel] = cfr_suite.cfr_optimize(model, onRxns, offRxns);
    %}

    arguments
        model (1,1) struct {mustHaveFields(model, {'S', 'rxns', 'lb', 'ub', 'c'})}
        onList {mustBeText}
        offList {mustBeText}
        options.listType (1,1) string {mustBeMember(options.listType, ["auto", "reactions", "genes"])} = "auto"
        options.onWeight {mustBeNumeric, mustBeNonnegative, mustBeLessThanOrEqual(options.onWeight, 10)} = 0.01
        options.onMinFlux {mustBeNumeric, mustBeNonnegative, mustBeLessThanOrEqual(options.onMinFlux, 1)} = 1e-3
        options.offWeight {mustBeNumeric, mustBeNonnegative, mustBeLessThanOrEqual(options.offWeight, 10)} = 0.01
        options.offTolerance {mustBeNumeric, mustBeNonnegative, mustBeLessThanOrEqual(options.offTolerance, 1)} = 0
        options.pfba (1,1) logical = true
        options.pfbaWeight (1,1) {mustBeNumeric, mustBeNonnegative, mustBeLessThanOrEqual(options.pfbaWeight, 10)} = 1e-6
        options.pfbaTolerance (1,1) {mustBeNumeric, mustBeNonnegative, mustBeLessThanOrEqual(options.pfbaTolerance, 1)} = 0
        options.gurobiParams (1,1) struct = struct()
        options.verbose (1,1) logical = true
        options.checkGurobiParams (1,1) logical = true
    end

%% RESOLVE GENE/REACTION LISTS AND SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%

    % Convert lists to reactions; keep per-reaction parameters aligned
    [onRxns, onParams] = resolve_list(model, onList, 'onList', options.listType, ...
        {options.onWeight, options.onMinFlux}, {'onWeight', 'onMinFlux'}, options.verbose);
    [offRxns, offParams] = resolve_list(model, offList, 'offList', options.listType, ...
        {options.offWeight, options.offTolerance}, {'offWeight', 'offTolerance'}, options.verbose);
    [onWeight, onMinFlux] = onParams{:};
    [offWeight, offTolerance] = offParams{:};

    % Consolidate contradicting reactions ("on" takes precedence)
    keep = ~ismember(offRxns, onRxns);
    offRxns = offRxns(keep);
    offWeight = offWeight(keep);
    offTolerance = offTolerance(keep);

    % Sort and remove duplicates, carrying parameters along
    [onRxns, ia] = unique(onRxns);
    onWeight = onWeight(ia);
    onMinFlux = onMinFlux(ia);
    [offRxns, ia] = unique(offRxns);
    offWeight = offWeight(ia);
    offTolerance = offTolerance(ia);

    % Gurobi parameters: drop invalid ones (as cfr_suite.py); don't display output
    gurobiParams = options.gurobiParams;
    if options.checkGurobiParams
        gurobiParams = cfr_suite.check_gurobi_params(gurobiParams, options.verbose);
    end
    gurobiParams.OutputFlag = 0;
    if ~any(strcmpi(fieldnames(gurobiParams), 'IntFeasTol'))
        gurobiParams.IntFeasTol = 1e-9;   % default: strict integrality for the big-M "on" constraints
    end

    % pFBA reactions: all reactions that are not "off"
    if options.pfba
        pfbaRxns = setdiff(model.rxns, offRxns);
        pfbaRxns = pfbaRxns(:);
    else
        pfbaRxns = cell(0, 1);
    end

    % No CFR constraints: the model is solved as plain FBA (as cfr_suite.py)
    if options.verbose && isempty(onRxns) && isempty(offRxns) && ~options.pfba
        fprintf('No CFR constraints detected: returning default solution\n');
    end

%% CONVERT TO GUROBI FORMAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S = sparse(double(model.S));   % gurobi requires a sparse double matrix
    [nMets, nRxns] = size(S);
    if isfield(model, 'b') && ~isempty(model.b)
        rhs = model.b(:);
    else
        rhs = zeros(nMets, 1);
    end
    if isfield(model, 'csense') && ~isempty(model.csense)
        sense = char(model.csense(:));
        sense(sense == 'E') = '=';
        sense(sense == 'L') = '<';
        sense(sense == 'G') = '>';
    else
        sense = repmat('=', nMets, 1);
    end
    bigM = 10000;

%% ADD CONSTRAINTS FOR "ON", "OFF" AND PFBA REACTIONS %%%%%%%%%%%%%%%%%%%%%

    % "on":          v - (eps1+M)*y1 >= -M    v + (eps1+M)*y2 <= M
    %                y1, y2 binary; maximize onWeight*(y1 + y2)
    % "off"/pFBA:    v + s >= -eps            v - r <= eps
    %                s, r >= 0; minimize weight*(s + r)
    % Rows/columns are interleaved per reaction to match cfr_suite.py:
    % [lower_r1; upper_r1; lower_r2; upper_r2; ...]
    slackRxns = [offRxns; pfbaRxns];
    slackWeight = [offWeight; repmat(options.pfbaWeight, numel(pfbaRxns), 1)];
    slackTolerance = [offTolerance; repmat(options.pfbaTolerance, numel(pfbaRxns), 1)];
    [~, onIdx] = ismember(onRxns, model.rxns);
    [~, slackIdx] = ismember(slackRxns, model.rxns);
    nOn = numel(onIdx);
    nSlack = numel(slackIdx);
    nAux = 2*nOn + 2*nSlack;

    % Interleave two column vectors: [a1; b1; a2; b2; ...]
    pair = @(a, b) reshape([a(:), b(:)].', [], 1);

    rows = (1:nAux)';
    Arxn = sparse(rows, [pair(onIdx, onIdx); pair(slackIdx, slackIdx)], 1, nAux, nRxns);
    Aaux = sparse(rows, rows, ...
        [pair(-(onMinFlux + bigM), onMinFlux + bigM); pair(ones(nSlack, 1), -ones(nSlack, 1))], nAux, nAux);

    gurobiModel = struct();
    gurobiModel.A = [S, sparse(nMets, nAux); Arxn, Aaux];
    gurobiModel.rhs = [rhs; pair(-bigM*ones(nOn, 1), bigM*ones(nOn, 1)); ...
        pair(-slackTolerance, slackTolerance)];
    gurobiModel.sense = [sense; repmat('><', 1, nOn + nSlack).'];
    gurobiModel.lb = [model.lb(:); zeros(nAux, 1)];
    gurobiModel.ub = [model.ub(:); ones(2*nOn, 1); 1000*ones(2*nSlack, 1)];
    gurobiModel.vtype = [repmat('C', nRxns, 1); repmat('B', 2*nOn, 1); repmat('C', 2*nSlack, 1)];
    gurobiModel.obj = [model.c(:); pair(onWeight, onWeight); -pair(slackWeight, slackWeight)];
    gurobiModel.modelsense = 'max';

%% DEFINE FUNCTION OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Solve and return a COBRA solution structure
    solverResult = gurobi(gurobiModel, gurobiParams);
    solution = struct();
    solution.solver = 'gurobi';
    solution.origStat = solverResult.status;
    switch solverResult.status
        case 'OPTIMAL'
            solution.stat = 1;
        case {'INFEASIBLE', 'INF_OR_UNBD'}
            solution.stat = 0;
        case 'UNBOUNDED'
            solution.stat = 2;
        case 'SUBOPTIMAL'
            solution.stat = 3;
        otherwise
            solution.stat = -1;
    end
    solution.onReactions = onRxns;
    solution.offReactions = offRxns;
    if solution.stat == 1
        solution.v = solverResult.x(1:nRxns);
        solution.f = model.c(:)' * solution.v;
        solution.cfrObjective = solverResult.objval;
        tol = 1e-9;
        solution.onSatisfied = abs(solution.v(onIdx)) >= onMinFlux - tol;
        solution.offSatisfied = abs(solution.v(slackIdx(1:numel(offRxns)))) <= offTolerance + tol;
    else
        warning("Gurobi status %s: returning NaN values.", solverResult.status)
        solution.v = NaN(nRxns, 1);
        solution.f = NaN;
        solution.cfrObjective = NaN;
        solution.onSatisfied = false(nOn, 1);
        solution.offSatisfied = false(numel(offRxns), 1);
    end
    solution.x = solution.v;
    solution = orderfields(solution, {'f', 'v', 'x', 'stat', 'origStat', 'solver', 'cfrObjective', ...
        'onReactions', 'offReactions', 'onSatisfied', 'offSatisfied'});

end

function [rxns, params] = resolve_list(model, list, label, listType, params, paramNames, verbose)
    % Convert a gene or reaction list to model reactions (as get_reactions
    % in cfr_suite.py) and expand parameters to one value per reaction
    ids = cellstr(list);
    ids = ids(:);
    if isempty(ids)
        rxns = cell(0, 1);
        params = cellfun(@(p) zeros(0, 1), params, 'UniformOutput', false);
        return
    end
    isRxn = ismember(ids, model.rxns);
    if isfield(model, 'genes')
        inGenes = ismember(ids, model.genes);
    else
        inGenes = false(size(ids));
    end
    isGene = inGenes & ~isRxn;
    if listType == "auto"
        if any(isRxn) && any(isGene)
            error('%s mixes %d reactions and %d genes: provide reaction IDs only or gene IDs only', ...
                label, sum(isRxn), sum(isGene));
        elseif any(isGene)
            listType = "genes";
        else
            listType = "reactions";
        end
    end
    if listType == "reactions"
        known = isRxn;
        kind = 'reactions';
    else
        known = inGenes;
        kind = 'genes';
    end
    if ~any(known)
        warning('No valid %s detected in %s', kind, label);
    elseif ~all(known)
        unknown = unique(ids(~known));
        warning('%d %s in %s are not in the model and are ignored (e.g. %s)', ...
            numel(unknown), kind, label, strjoin(unknown(1:min(5, end))', ', '));
    end

    if listType == "reactions"
        params = cellfun(@(p, name) expand_param(p, numel(ids), name), ...
            params, paramNames, 'UniformOutput', false);
        rxns = ids(known);
        params = cellfun(@(p) p(known), params, 'UniformOutput', false);
    else
        if any(cellfun(@numel, params) > 1)
            error('With gene lists, use scalar values for %s', strjoin(paramNames, ' and '));
        end
        if any(known)
            [~, ~, rxns, ~] = deleteModelGenes(model, unique(ids(known)));
            rxns = rxns(:);
            if verbose
                fprintf('Determined %s reactions from rxnGeneMat\n', label);
            end
        else
            rxns = cell(0, 1);
        end
        params = cellfun(@(p, name) expand_param(p, numel(rxns), name), ...
            params, paramNames, 'UniformOutput', false);
    end
end

function p = expand_param(p, n, name)
    % Expand a scalar to n values; check vector lengths
    if numel(p) == 1
        p = repmat(p, n, 1);
    elseif numel(p) ~= n
        error('%s must be a scalar or have one value per list entry (%d)', name, n);
    end
    p = double(p(:));
end

function mustHaveFields(s, f)
    % Validate that a model structure has the required fields
    missing = f(~isfield(s, f));
    if ~isempty(missing)
        error('model is missing required fields: %s', strjoin(missing, ', '));
    end
end
