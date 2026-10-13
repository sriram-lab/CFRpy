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
%
% MATLAB version compatibility
%   The package needs MATLAB R2020b or later (R2022b or later recommended).
%   Minimum release per function (own syntax requirements):
%     cfr_optimize        R2020b  arguments block with name-value options, mustBeText
%     apply_cfr           R2020b  calls cfr_optimize (own code uses older syntax)
%     check_gurobi_params R2019b  arguments block
%     summarize           R2019b  arguments block
%     get_fluxes          R2019b  arguments block with name-value options
%   Examples written as Name=Value (e.g. onWeight=1) need R2021a or later; in
%   R2020b use the 'name', value form (e.g. 'onWeight', 1).
%   No MATLAB toolboxes are required (base MATLAB only).
%
% Required external software (not included with MATLAB)
%   Gurobi (with MATLAB interface)   all functions that solve a model:
%                                    cfr_optimize, apply_cfr, check_gurobi_params
%   COBRA Toolbox                    cfr_optimize and apply_cfr with gene
%                                    lists (uses deleteModelGenes); not needed
%                                    when only reaction IDs are given to
%                                    cfr_optimize, nor for summarize/get_fluxes
%
%   Install Gurobi
%     1. Download and install Gurobi from https://www.gurobi.com (academic
%        licences are free; they must be renewed annually).
%     2. Activate the licence with the grbgetkey command from the Gurobi
%        download page.
%     3. In MATLAB, run the setup script of the Gurobi MATLAB interface:
%          >> cd('<gurobi_home>/<platform>/matlab')   % e.g. C:\gurobi1303\win64\matlab
%          >> gurobi_setup
%          >> savepath                                % optional: keep on the path
%     4. Check that it is found: >> which gurobi
%     Use a Gurobi release that supports your MATLAB release (see the Gurobi
%     documentation), and the same major Gurobi version as in Python when
%     comparing results.
%
%   Install the COBRA Toolbox
%     1. Download it (git is the easiest way; or download the ZIP from
%        https://github.com/opencobra/cobratoolbox):
%          git clone --depth=1 https://github.com/opencobra/cobratoolbox.git
%     2. In MATLAB, initialize it (updateToolbox = false):
%          >> cd('<path>/cobratoolbox')
%          >> initCobraToolbox(false)
%     3. Check that it is found: >> which deleteModelGenes
%     initCobraToolbox must be run in each new MATLAB session (or call it
%     from startup.m).
