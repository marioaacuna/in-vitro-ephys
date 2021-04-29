function run_in_vitro_ephys_analysis()
% Set global variables
global GC
% Read general_configs
GC = general_configs();
%%
addpath(genpath('C:\Repositories\in-vitro-ephys'))
AP_launcher_app()
