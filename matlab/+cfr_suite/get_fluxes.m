function fluxes = get_fluxes(results, model, options)
    % GET_FLUXES Extract flux solutions of apply_cfr results into a table.
    %
    % DESCRIPTION
    % Collects the flux solution of each condition into a reactions x
    % conditions table, optionally preceded by reaction metadata (reaction
    % names and subsystems). Equivalent to get_fluxes in cfr_suite.py.
    %
    % Author:   Carolina H. Chung
    % Created:  October 8, 2026

    % I/O
    %{
    INPUTS:
        1. results:     output of cfr_suite.apply_cfr
           -> type:     structure array (condition, onList, offList,
                        solution)

        2. model:       the model used in apply_cfr; its rxns are used as
                        row names (and rxnNames/subSystems as metadata).
                        Without a model, rows are in reaction order with
                        no names (optional)
           -> type:     MATLAB structure

    OPTIONS (name-value pairs):
        metadata:       add Reaction (model.rxnNames) and Subsystem
                        (model.subSystems) before the flux columns;
                        requires model
           -> default:  false

    OUTPUTS:
        1. fluxes:      one row per reaction, one variable per condition
                        (preceded by Reaction and Subsystem with metadata)
           -> type:     table

    EXAMPLE USAGE:
        >> results = cfr_suite.apply_cfr(model, data);
        >> fluxes = cfr_suite.get_fluxes(results, model);
        >> fluxes = cfr_suite.get_fluxes(results, model, metadata=true);
    %}

    arguments
        results struct
        model struct = struct()
        options.metadata (1,1) logical = false
    end

    % Collect flux solutions (reactions x conditions)
    results = results(:);
    v = arrayfun(@(r) r.solution.v(:), results, 'UniformOutput', false);
    nRxns = cellfun(@numel, v);
    if any(nRxns ~= nRxns(1))
        error('Flux solutions have different lengths: results come from different models');
    end
    hasModel = isfield(model, 'rxns');
    if hasModel && numel(model.rxns) ~= nRxns(1)
        error('model has %d reactions, but flux solutions have %d', numel(model.rxns), nRxns(1));
    end

    % Create flux table
    fluxes = array2table([v{:}], 'VariableNames', {results.condition});
    if hasModel
        fluxes.Properties.RowNames = model.rxns(:);
    end

    % Add metadata (if prompted)
    if options.metadata
        if ~hasModel
            error('Provide the model to add reaction metadata');
        end
        n = numel(model.rxns);
        reaction = repmat({''}, n, 1);
        if isfield(model, 'rxnNames')
            reaction = cellstr(model.rxnNames(:));
        end
        subsystem = repmat({''}, n, 1);
        if isfield(model, 'subSystems')
            subsystem = cellfun(@subsystem_name, model.subSystems(:), 'UniformOutput', false);
        end
        fluxes = [table(reaction, subsystem, 'VariableNames', {'Reaction', 'Subsystem'}), fluxes];
    end

end

function name = subsystem_name(s)
    % Subsystem as text (COBRA models store one or more subsystems per reaction)
    if isempty(s)
        name = '';
    else
        name = strjoin(cellstr(s), '; ');
    end
end
