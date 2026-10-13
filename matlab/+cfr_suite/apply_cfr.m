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

        progressBar:    show a waitbar with the progress across conditions
                        and an estimate of the time remaining (skipped when
                        MATLAB has no display, e.g. with -nodisplay)
           -> type:     logical (or 0/1)
           -> default:  true

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

        3. Without the progress bar:
           >> results = cfr_suite.apply_cfr(model, data, progressBar=false);
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
    progressBar = true;
    isBar = strcmpi(names, 'progressBar');
    if any(isBar)
        progressBar = values{find(isBar, 1, 'last')};
    end
    names = names(~isBar);
    values = values(~isBar);
    if ~(islogical(progressBar) || isnumeric(progressBar)) || ~isscalar(progressBar) ...
            || ~ismember(progressBar, [0 1])
        error('Provide true or false for progressBar');
    end
    progressBar = logical(progressBar);

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

    % Show progress with a waitbar (needs a display); it closes when
    % finished or if an error interrupts the loop
    conditions = data.Properties.VariableNames;
    nConditions = numel(conditions);
    bar = [];
    if progressBar && usejava('awt') && nConditions > 0
        bar = waitbar(0, sprintf('Applying CFR across %d conditions (0%%)', nConditions), 'Name', 'CFR');
        closeProgress = onCleanup(@() close_progress(bar));
    end

    % Determine solutions
    results = struct('condition', conditions, 'onList', [], 'offList', [], 'solution', []);
    startTime = tic;
    for k = 1:numel(conditions)
        levels = data.(conditions{k});
        results(k).onList = genes(levels > thresh(2));
        results(k).offList = genes(levels < thresh(1));
        results(k).solution = cfr_suite.cfr_optimize(model, results(k).onList, results(k).offList, cfrArgs{:});
        if isgraphics(bar)   % false if there is no bar or the user closed it
            remaining = toc(startTime) / k * (nConditions - k);
            waitbar(k/nConditions, bar, ...
                sprintf('Applying CFR across %d conditions (%.0f%%)\n%s remaining', ...
                nConditions, 100*k/nConditions, format_duration(remaining)));
        end
    end

end

function text = format_duration(seconds)
    % Format a duration in seconds as hours:minutes:seconds (H:MM:SS)
    seconds = round(seconds);
    text = sprintf('%d:%02d:%02d', floor(seconds/3600), floor(mod(seconds, 3600)/60), mod(seconds, 60));
end

function close_progress(bar)
    % Close the waitbar if it is still open
    if isgraphics(bar)
        close(bar);
    end
end
