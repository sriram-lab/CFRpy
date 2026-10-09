function summary = summarize(results)
    % SUMMARIZE Summarize apply_cfr results into a table.
    %
    % DESCRIPTION
    % Collects the number of up-regulated genes (N_up), down-regulated genes
    % (N_down) and the objective value (Objective) of each condition.
    % Equivalent to summarize in cfr_suite.py.
    %
    % Author:   Carolina H. Chung
    % Created:  October 8, 2026

    % I/O
    %{
    INPUTS:
        1. results:     output of cfr_suite.apply_cfr
           -> type:     structure array (condition, onList, offList,
                        solution)

    OUTPUTS:
        1. summary:     one row per condition (row names), with variables
                        N_up, N_down and Objective (solution.f)
           -> type:     table

    EXAMPLE USAGE:
        >> results = cfr_suite.apply_cfr(model, data);
        >> summary = cfr_suite.summarize(results)
    %}

    arguments
        results struct
    end

    % Collect gene counts and objective values
    results = results(:);
    nUp = arrayfun(@(r) numel(r.onList), results);
    nDown = arrayfun(@(r) numel(r.offList), results);
    objective = arrayfun(@(r) r.solution.f, results);

    % Create summary table
    summary = table(nUp, nDown, objective, 'VariableNames', {'N_up', 'N_down', 'Objective'}, ...
        'RowNames', {results.condition}');

end
