function run_in_vitro_ephys_analysis()
% Set global variables
global GC
% Read general_configs
GC = general_configs();
%%
addpath(genpath(cd()))
AP_launcher_app()
