function [fluxstate, grate, solverobj, solg1, model] = cfr_optimize_legacy(...
    model1, onreactions, offreactions, ...
    kappa, rho, epsilon, mode, epsilon2, minfluxflag, kappa2, params)

    % DESCRIPTION
    % This function simulates flux through reactions annotated in a
    % genome-scale metabolic model (GEM) by constraining the model using
    % information on "on" and "off" reactions (can be based directly or by
    % differentially expressed gene (DEG) data.
    %
    % Based on constrain_flux_regulation.m, with these fixes:
    % - per-reaction parameters stay paired with their reactions after
    %   contradicting reactions are removed and the lists are sorted
    %   (previously, leftover "off" weights were assigned to pFBA reactions)
    % - model1.csense is applied (previously always ignored)
    % - reaction lists of any orientation are accepted
    % - reactions not in the model are ignored with a warning (previously,
    %   they were silently ignored)
    % - a full (non-sparse) model1.S is accepted
    % - params.IntFeasTol defaults to 1e-9, so "on" binaries cannot be
    %   credited without flux under the big-M constraints
    % - non-optimal Gurobi solutions return NaN values
    % - "on" and "off" constraints are built in one step (faster)
    %
    % STEPS
    % 1. Inspect inputs and set default values
    % 2. Convert model to gurobi format
    % 3. Add constraints for "on" and "off" reactions
    % 4. Define function outputs
    %
    % Author:   Sriram Chandrasekaran
    % Created:  October 12, 2010
    % Updates:  October 6, 2026 (Carolina H. Chung)
    %           July 14, 2026 (Carolina H. Chung)
    %           September 20, 2020 (Carolina H. Chung)
    %           November 6, 2019 (Scott Campit)
    %           November 8, 2018 (Fangzhou Shen)

    % I/O
    %{
    REQUIRED INPUTS:
        1. model1:          genome-scale metabolic model (GEM)
           -> type:         MATLAB structure

        2. onreactions:     list of "on" reactions (when 'mode' = 1) or
                            list of up-regulated genes (when 'mode' = 0)
           -> type:         cell array
           -> dim:          m x 1 (m = # of "on" reactions/up-reg genes)

        3. offreactions:    list of  "off" reactions (when 'mode' = 1) or
                            list of down-regulated genes (when 'mode' = 0)
           -> type:         cell array
           -> dim:          m x 1 (m = # of "off" reactions/down-reg genes)

    OPTIONAL INPUTS (in argument order):
        1. kappa:       relative weight for "off" reactions
           -> default:  kappa = 1
           -> note:     use small value (e.g. 1e-3) for large networks

        2. rho:         relative weight for "on" reactions
           -> default:  rho = 1
           -> note:     use small value (e.g. 1e-3) for large networks

        3. epsilon:     minimum flux for "on" reactions
           -> default:  epsilon = 1e-3

        4. mode:        specifies whether inputs are reaction or gene lists
                        (mode = 0 for genes     mode = 1 for reactions)
           -> default:  mode = 1

        5. epsilon2:    flux tolerance for "off" reactions
           -> default:  epsilon2 = 0

        6. minfluxflag: Boolean to minimize flux sum through all reactions
                        (i.e. apply PFBA)
           -> default:  minfluxflag = true

        7. kappa2:      relative weight for pFBA reactions
           -> default:  kappa2 = 1e-6

        8. params:      parameter values to pass to gurobi
           -> default:  params.OutputFlag = 0, params.IntFeasTol = 1e-9

        kappa, rho, epsilon and epsilon2 are scalars or vectors with one
        value per input reaction (mode = 1). With gene lists (mode = 0),
        use scalars.

    OUTPUTS:
        1. fluxstate:       flux solution through reactions in GEM
        2. grate:           simulated growth rate
        3. solverobj:       flux through GEM objective function
        4. solg1:           gurobi solution structure
        5. model:           gurobi model with CFR constraints

    EXAMPLE USAGE:
        1. Default usage:
           >> [fluxstate, grate, solverobj] = cfr_optimize_legacy(...
                    model1, onreaction_list, offreaction_list)

        2. Use gene lists instead with defaults:
           >> [fluxstate, grate, solverobj] = cfr_optimize_legacy(...
                    model1, up_genes, down_genes, [], [], [], 0)
    %}

%% INSPECT INPUTS AND SET DEFAULT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % mode (reactions or genes)
    if (~exist('mode','var')) || (isempty(mode))
        mode = 1;
    end
    if mode == 0
        [~, ~, onreactions, ~] =  deleteModelGenes(model1, onreactions);
        [~, ~, offreactions, ~] =  deleteModelGenes(model1, offreactions);
    end
    onreactions = onreactions(:);
    offreactions = offreactions(:);

    % Function parameters (scalars, or one value per reaction)
    if (~exist('epsilon','var')) || (isempty(epsilon))
        epsilon = 1E-3;
    end
    if (~exist('rho','var')) || (isempty(rho))
        rho = 1;
    end
    if (~exist('kappa','var')) || (isempty(kappa))
        kappa = 1;
    end
    if (~exist('epsilon2','var')) || (isempty(epsilon2))
        epsilon2 = 0;
    end
    if (~exist('minfluxflag','var')) || (isempty(minfluxflag))
        minfluxflag = true;
    else
        minfluxflag = logical(minfluxflag);
    end
    if (~exist('kappa2','var')) || (isempty(kappa2))
        kappa2 = 1E-6;
    end
    epsilon = expand_param(epsilon, onreactions, 'epsilon');
    rho = expand_param(rho, onreactions, 'rho');
    kappa = expand_param(kappa, offreactions, 'kappa');
    epsilon2 = expand_param(epsilon2, offreactions, 'epsilon2');

    % Drop reactions not in the model
    unknown = setdiff([onreactions; offreactions], model1.rxns);
    if ~isempty(unknown)
        warning('%d reactions are not in model1.rxns and are ignored (e.g. %s)', ...
            numel(unknown), strjoin(unknown(1:min(5, end))', ', '));
    end
    keep = ismember(onreactions, model1.rxns);
    onreactions = onreactions(keep);
    epsilon = epsilon(keep);
    rho = rho(keep);

    % Consolidate contradicting reactions ("on" takes precedence)
    keep = ismember(offreactions, model1.rxns) & ...
        ~ismember(offreactions, onreactions);
    offreactions = offreactions(keep);
    kappa = kappa(keep);
    epsilon2 = epsilon2(keep);

    % Sort and remove duplicates, carrying parameters along
    [onreactions, ia] = unique(onreactions);
    epsilon = epsilon(ia);
    rho = rho(ia);
    [offreactions, ia] = unique(offreactions);
    kappa = kappa(ia);
    epsilon2 = epsilon2(ia);

    % Gurobi parameters: don't display output; strict integrality for big-M
    if (~exist('params','var')) || (isempty(params))
        params = struct();
    end
    if ~isfield(params, 'OutputFlag')
        params.OutputFlag = 0;
    end
    if ~isfield(params, 'IntFeasTol')
        params.IntFeasTol = 1e-9;
    end

    % Apply PFBA (if prompted)
    if minfluxflag
        pfbareactions = setdiff(model1.rxns, offreactions);
        kappa = [kappa; repmat(kappa2, numel(pfbareactions), 1)];
        epsilon2 = [epsilon2; zeros(numel(pfbareactions), 1)];
        offreactions = [offreactions; pfbareactions(:)];
    end

%% CONVERT TO GUROBI FORMAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Re-define model
    model = model1;
    model.A = sparse(double(model1.S));   % gurobi requires a sparse double matrix
    model.obj = model1.c;
    if isfield(model1, 'b') && ~isempty(model1.b)
        model.rhs = model1.b(:);
    else
        model.rhs = zeros(size(model1.S, 1), 1);
    end
    if isfield(model1, 'csense') && ~isempty(model1.csense)
        model.sense = char(model1.csense(:));
        model.sense(model.sense == 'E') = '=';
        model.sense(model.sense == 'L') = '<';
        model.sense(model.sense == 'G') = '>';
    else
        model.sense = repmat('=', [size(model1.S, 1), 1]);
    end
    model.lb = model1.lb;
    model.ub = model1.ub;
    model.vtype = repmat('C', size(model1.S, 2), 1);
    model.modelsense = 'max';
    M = 10000;
    nrxns = length(model1.rxns);

%% ADD CONSTRAINTS FOR "ON" AND "OFF" REACTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % "on":          v - (eps+M)*y1 >= -M    v + (eps+M)*y2 <= M
    %                y1, y2 binary; maximize rho*(y1 + y2)
    % "off"/pFBA:    v + s >= -eps2          v - r <= eps2
    %                s, r >= 0; minimize kappa*(s + r)
    [~, onidx] = ismember(onreactions, model1.rxns);
    [~, offidx] = ismember(offreactions, model1.rxns);
    non = numel(onidx);
    noff = numel(offidx);
    nnew = 2*non + 2*noff;
    [m, n] = size(model1.S);

    rows = (1:nnew)';
    Arxn = sparse(rows, [onidx; onidx; offidx; offidx], 1, nnew, n);
    Aaux = sparse(rows, rows, ...
        [-(epsilon + M); (epsilon + M); ones(noff, 1); -ones(noff, 1)], nnew, nnew);

    model.A = [model.A, sparse(m, nnew); Arxn, Aaux];
    model.rhs = [model.rhs; -M*ones(non, 1); M*ones(non, 1); -epsilon2; epsilon2];
    model.sense = [model.sense; repmat('>', non, 1); repmat('<', non, 1); ...
        repmat('>', noff, 1); repmat('<', noff, 1)];
    model.lb = [model1.lb(:); zeros(nnew, 1)];
    model.ub = [model1.ub(:); ones(2*non, 1); 1000*ones(2*noff, 1)];
    model.vtype = [repmat('C', n, 1); repmat('B', 2*non, 1); repmat('C', 2*noff, 1)];
    model.obj = [model1.c(:); rho; rho; -kappa; -kappa];

%% DEFINE FUNCTION OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    solg1 = gurobi(model, params);
    if strcmp(solg1.status, 'OPTIMAL')
        fluxstate = solg1.x(1:nrxns);
        grate = dot(fluxstate, model1.c);
        solverobj = solg1.objval;
    else
        warning("Gurobi status %s: returning NaN values.", solg1.status)
        fluxstate = NaN .* ones(nrxns, 1);
        grate = NaN;
        solverobj = NaN;
    end

end

function p = expand_param(p, rxns, name)
    % Expand a scalar to one value per reaction; check vector lengths
    if numel(p) == 1
        p = repmat(p, numel(rxns), 1);
    elseif numel(p) ~= numel(rxns)
        error('%s must be a scalar or have one value per reaction (%d)', name, numel(rxns));
    end
    p = p(:);
end
