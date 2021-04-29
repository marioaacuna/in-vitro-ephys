function run_in_vitro_ephys_analysis()
% Set global variables
global GC
% Read general_configs
GC = general_configs();
%%
addpath(genpath('M:\Mario\Fede\Ephys_in_vitro_analysis'))
AP_launcher_app()
