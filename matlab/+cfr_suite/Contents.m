% CFR_SUITE MATLAB implementation of constrain flux regulation (CFR)
%
% Equivalent of the Python module cfr_suite.py. Call each function with the
% package prefix, e.g. cfr_suite.apply_cfr(model, data).
%
% Functions
%   cfr_optimize        - Optimize a COBRA model via CFR
%   apply_cfr           - Apply CFR across a set of conditions
%   summarize           - Summarize apply_cfr results (N_up, N_down, Objective)
%   get_fluxes          - Flux table (reactions x conditions) of apply_cfr results
%   check_gurobi_params - Validate a structure of Gurobi parameters
%
% Legacy functions (matlab/legacy, outside the package)
%   cfr_optimize_legacy       - constrain_flux_regulation with bug fixes
%   constrain_flux_regulation - original CFR implementation
%
% Setup
%   addpath('<CFRpy>/matlab')            % cfr_suite package
%   addpath('<CFRpy>/matlab/legacy')     % legacy functions (optional)
