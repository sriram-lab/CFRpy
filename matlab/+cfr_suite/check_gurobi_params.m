function params = check_gurobi_params(gurobiParams, verbose)
    % CHECK_GUROBI_PARAMS Validate a structure of Gurobi parameters.
    %
    % DESCRIPTION
    % Each parameter is checked against the parameter names of the Gurobi
    % installation (gurobi ignores unknown names without notice) and tested
    % on a trivial model (gurobi rejects invalid values with an error).
    % Unknown names, rejected values and duplicated names (ignoring case)
    % raise a warning; invalid parameters are dropped and the last value of
    % a duplicated name is kept. Valid parameters are printed when verbose
    % is true. Equivalent to check_gurobi_params in cfr_suite.py.
    %
    % Author:   Carolina H. Chung
    % Created:  October 8, 2026

    % I/O
    %{
    INPUTS:
        1. gurobiParams:    gurobi parameters (field names are parameter
                            names, in any case)
           -> type:         MATLAB structure
           -> default:      struct()

        2. verbose:         print the name and value of each valid parameter
           -> type:         logical
           -> default:      true

    OUTPUTS:
        1. params:          valid parameters, with gurobi parameter names
           -> type:         MATLAB structure

    WARNINGS:
        cfr_optimize:invalidGurobiParam     unknown name or rejected value
        cfr_optimize:duplicateGurobiParam   name given more than once

    EXAMPLE USAGE:
        >> params = cfr_suite.check_gurobi_params(struct('intfeastol', 1e-9, 'NotAParam', 1));
        Warning: Ignoring invalid Gurobi parameter NotAParam: Unknown parameter 'NotAParam'
        Set Gurobi parameter IntFeasTol to 1e-09
    %}

    arguments
        gurobiParams (1,1) struct = struct()
        verbose (1,1) logical = true
    end

    % Validate parameter names and values
    validNames = gurobi_param_names();
    probe = struct('A', sparse(1, 1), 'rhs', 0, 'sense', '<', 'obj', 0, 'lb', 0, 'ub', 0);
    params = struct();
    given = struct();
    keys = fieldnames(gurobiParams);
    for k = 1:numel(keys)
        key = keys{k};
        value = gurobiParams.(key);
        if isempty(validNames)
            name = key;   % parameter list not found: check values only
        else
            name = validNames(strcmpi(key, validNames));
            if isempty(name)
                warning('cfr_optimize:invalidGurobiParam', ...
                    'Ignoring invalid Gurobi parameter %s: Unknown parameter ''%s''', key, key);
                continue
            end
            name = name{1};
        end
        p = struct();
        p.(name) = value;
        p.OutputFlag = 0;
        if strcmp(name, 'OutputFlag')
            p.OutputFlag = value;
        end
        try
            evalc('gurobi(probe, p);');   % evalc hides any log the probe prints
        catch err
            reason = regexprep(strtrim(err.message), '^Gurobi error \d+: ', '');
            warning('cfr_optimize:invalidGurobiParam', ...
                'Ignoring invalid Gurobi parameter %s: %s', key, reason);
            continue
        end
        if isfield(params, name)
            warning('cfr_optimize:duplicateGurobiParam', ...
                'Duplicated Gurobi parameter %s (given as %s and %s): using %s', name, given.(name), key, key);
        end
        params.(name) = value;
        given.(name) = key;
    end

    % Print valid parameters (if prompted)
    if verbose
        names = fieldnames(params);
        for k = 1:numel(names)
            fprintf('Set Gurobi parameter %s to %s\n', names{k}, format_value(params.(names{k})));
        end
    end

end

function names = gurobi_param_names()
    % Gurobi parameter names, read once from gurobi_c.h of the installation
    % that provides the gurobi MEX file (empty if not found)
    persistent cache
    if isempty(cache)
        cache = {{}};
        root = fileparts(fileparts(which('gurobi')));
        header = fullfile(root, 'include', 'gurobi_c.h');
        if isfile(header)
            tok = regexp(fileread(header), '#define GRB_(?:INT|DBL|STR)_PAR_\w+\s+"(\w+)"', 'tokens');
            cache = {cellfun(@(t) t{1}, tok, 'UniformOutput', false)'};
        end
    end
    names = cache{1};
end

function s = format_value(value)
    % Format a parameter value as in cfr_suite.py ('%g' for numbers)
    if (isnumeric(value) || islogical(value)) && isscalar(value)
        s = sprintf('%g', value);
    elseif isnumeric(value) || islogical(value)
        s = mat2str(value);
    else
        s = char(string(value));
    end
end
