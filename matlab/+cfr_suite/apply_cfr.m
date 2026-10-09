function results = apply_cfr(model, data, varargin)
    % APPLY_CFR Apply CFR across a set of conditions.
    %
    % DESCRIPTION
    % For each condition (column) of an omics data table, genes above the
    % upper threshold are used as "on" genes and genes below the lower
    % threshold as "off" genes in cfr_optimize. Equivalent to apply_cfr in
    % cfr_suite.py.
    %
    % Author:   Carolina H. Chung
    % Created:  October 8, 2026

    % I/O
    %{
    REQUIRED INPUTS:
        1. model:       genome-scale metabolic model (GEM) with a genes
                        field (see cfr_optimize for other fields)
           -> type:     MATLAB structure

        2. data:        omics data used to infer "on" (up-regulated) and
                        "off" (down-regulated) genes; one numeric variable
                        per condition, row names matching model.genes
           -> type:     table

    OPTIONAL INPUTS:
        3. thresh:      [lower, upper] thresholds for "off" (< lower) and
                        "on" (> upper) genes, given as the third argument
                        or as the name-value pair 'thresh'
           -> default:  [-2, 2]

        Any other name-value pairs are passed to cfr_optimize (e.g.
        onWeight, offWeight, pfbaWeight, gurobiParams, verbose).
        gurobiParams is validated once for all conditions.

        progressBar:    progress across conditions, shown with progressbar
                        (MATLAB File Exchange) when it is on the path:
                        1 = own bar, n > 1 = update bar n of a progressbar
                        opened by the caller (e.g. an outer loop), 0 or
                        false = no bar
           -> default:  1

    OUTPUTS:
        1. results:     one element per condition (order of data columns)
           -> type:     structure array with fields
              condition:    condition name (data column)
              onList:       up-regulated genes
              offList:      down-regulated genes
              solution:     cfr_optimize solution

    EXAMPLE USAGE:
        1. Default thresholds:
           >> results = cfr_suite.apply_cfr(model, data);
           >> cfr_suite.summarize(results)

        2. Custom thresholds and CFR options:
           >> results = cfr_suite.apply_cfr(model, data, [-1.5 1.5], onWeight=1, ...
                    offWeight=1, pfbaWeight=1e-4, verbose=false);

        3. Sub-bar inside an outer loop (the outer loop creates both bars):
           >> progressbar('Objectives', 'Conditions');
           >> for i = 1:numel(models)
           >>     results{i} = cfr_suite.apply_cfr(models{i}, data, progressBar=2);
           >>     progressbar(i/numel(models));
           >> end
    %}

    % Check model
    if ~isstruct(model) || ~isfield(model, 'genes')
        error('Invalid COBRA model: requires a genes field');
    end

    % Check data
    if ~istable(data) || isempty(data.Properties.RowNames)
        error('Provide data as a table with gene IDs as row names');
    end

    % Separate thresh from the name-value pairs passed to cfr_optimize
    thresh = [-2 2];
    if ~isempty(varargin) && isnumeric(varargin{1})
        thresh = varargin{1};
        varargin(1) = [];
    end
    if mod(numel(varargin), 2) ~= 0
        error('Provide name-value pairs after thresh');
    end
    names = cellfun(@char, varargin(1:2:end), 'UniformOutput', false);
    values = varargin(2:2:end);
    isThresh = strcmpi(names, 'thresh');
    if any(isThresh)
        thresh = values{find(isThresh, 1, 'last')};
    end
    names = names(~isThresh);
    values = values(~isThresh);

    % Separate progressBar (apply_cfr option, not passed to cfr_optimize)
    progressBar = 1;
    isBar = strcmpi(names, 'progressBar');
    if any(isBar)
        progressBar = double(values{find(isBar, 1, 'last')});
    end
    names = names(~isBar);
    values = values(~isBar);
    if ~isscalar(progressBar) || progressBar < 0 || progressBar ~= round(progressBar)
        error('Provide a non-negative integer (or logical) for progressBar');
    end

    % Check thresh
    if ~isnumeric(thresh) || numel(thresh) ~= 2
        error('Provide a numeric vector of two values for thresh');
    end

    % Validate gurobi parameters once for all conditions
    isVerbose = strcmpi(names, 'verbose');
    verbose = true;
    if any(isVerbose)
        verbose = values{find(isVerbose, 1, 'last')};
    end
    isParams = strcmpi(names, 'gurobiParams');
    gurobiParams = struct();
    if any(isParams)
        gurobiParams = values{find(isParams, 1, 'last')};
    end
    gurobiParams = cfr_suite.check_gurobi_params(gurobiParams, verbose);
    keep = ~isParams & ~strcmpi(names, 'checkGurobiParams');
    names = [names(keep), {'gurobiParams', 'checkGurobiParams'}];
    values = [values(keep), {gurobiParams, false}];
    cfrArgs = reshape([names; values], 1, []);

    % Filter data
    genes = data.Properties.RowNames;
    if ~any(ismember(genes, model.genes))
        warning('Index of input data does not match COBRA model gene names');
    end
    data = data(ismember(genes, model.genes), :);
    genes = data.Properties.RowNames;

    % Show progress (if progressbar from MATLAB File Exchange is installed);
    % an own bar closes when finished or if an error interrupts the loop,
    % a sub-bar (progressBar > 1) leaves the caller's progressbar open
    conditions = data.Properties.VariableNames;
    showProgress = progressBar > 0 && exist('progressbar', 'file') == 2 && ~isempty(conditions);
    if showProgress && progressBar == 1
        progressbar(sprintf('Applying CFR across %d conditions', numel(conditions)));
        closeProgress = onCleanup(@() progressbar(1));
    elseif showProgress
        update_progress(progressBar, 0);
    end

    % Determine solutions
    results = struct('condition', conditions, 'onList', [], 'offList', [], 'solution', []);
    for k = 1:numel(conditions)
        levels = data.(conditions{k});
        results(k).onList = genes(levels > thresh(2));
        results(k).offList = genes(levels < thresh(1));
        results(k).solution = cfr_suite.cfr_optimize(model, results(k).onList, results(k).offList, cfrArgs{:});
        if showProgress
            update_progress(progressBar, k/numel(conditions));
        end
    end

end

function update_progress(bar, fraction)
    % Update one bar of progressbar, leaving the other bars unchanged
    args = repmat({[]}, 1, bar);
    args{bar} = fraction;
    progressbar(args{:});
end
