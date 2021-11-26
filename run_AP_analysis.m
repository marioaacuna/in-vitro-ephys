function DONE =  run_AP_analysis(experimenter_ID, recording_date, do_plotting, varargin)
clc
DONE = 0;
p = varargin;
overwrite = p{7};
% Set global variables
global GC
% Read general_configs
GC = general_configs();
%%
GC.experimenter_ID = experimenter_ID;
data_path_root =  GC.raw_data_root_path.(experimenter_ID); % This needs to be added on GC
toolboxes_to_use = {'Igor2Matlab'};
toggle_toolbox(toolboxes_to_use, 'on')
switch experimenter_ID
    case {'Federica', 'Fede_setup_1', 'Fede_setup_1_5HT7','Niels'}
        data_path = os.path.join(data_path_root, recording_date{1}, 'traces');
    case 'Liselot'
        data_path = os.path.join(data_path_root, '2. opto-5HT7', recording_date{1}, [recording_date{1},'.data']);
    case 'Sri'
        data_path = os.path.join(data_path_root, recording_date{1});
    otherwise
        keyboard 
end
% data_path = 'M:\Mario\Fede\traces'; % To be changed later
data_dir = dir(data_path);
% Pick animals
ending = '.ibw';
names ={data_dir.name};
files_in_folder = names(endsWith(names, ending));
files_in_folder = natsort(files_in_folder);

% Isolate the files to analyze
str_exptr = GC.string_file_selection.(experimenter_ID);
if any(strcmp(experimenter_ID, {'Niels', 'Sri', 'Fede_setup_1'})) 
   is_Amp = cell2mat(cellfun(@(x) sum(ismember(x,str_exptr)) == length(str_exptr)+1 && ~endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
else
    is_Amp = cell2mat(cellfun(@(x) sum(ismember(x,str_exptr)) == length(str_exptr) && ~endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
end
files_to_take = files_in_folder(is_Amp);
% is_outwave = cell2mat(cellfun(@(x) sum(ismember(x,str_exptr)) == length(str_exptr) && endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
is_outwave = cell2mat(cellfun(@(x) endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
outwaves_files = files_in_folder(is_outwave); 
outwave_identifier = GC.set_string_outwave_selection.(experimenter_ID);
%% Run analysis
FR_AP = [];%NaN(length(GC.current_steps.(experimenter_ID)),1);
amp_AP = [];
width_AP = [];
rest_Vm = [];
input_R = [];
Firing_threshold = [];
SAG_r = [];
SAG_d = [];
A_bump = [];
Vel_depo = [];
Vel_repo = [];
Tau_mb = [];
Area_phase = [];
AP_phase = [];
DER_phase =[];
Adapt_index = [];
Burst_freq = [];
% do_plotting = 0;
names = cell(length(files_to_take),1);
for i_exp = 1:length(files_to_take)
    try
        this_exp = files_to_take{i_exp};
        file_to_read = os.path.join(data_path,this_exp);
        D = IBWread(file_to_read); % read Voltage traces
        str_V = strsplit(this_exp,str_exptr);% first characters of the voltage filethat need to match to the Current file
        start_str_v = [str_V{1}, outwave_identifier] ; 
        outwave_file_start= char(outwaves_files(startsWith(outwaves_files, start_str_v)));
        % check for correct identifier
        if isempty(outwave_file_start)
            start_str_v = [str_V{1}, str_exptr] ; 
            outwave_file_start= char(outwaves_files(startsWith(outwaves_files, start_str_v)));       
        end
        outwave_file = os.path.join(data_path,outwave_file_start);
        O =  IBWread(outwave_file); % read Current traces
        % take only data that has current steps
        if size(O.y, 2) ~= size(D.y, 2)
%             disp([this_exp, ' did not correspond to I/O curve'])
            to_remove = size(O.y, 2) - size(D.y, 2);
            I_traces = O.y(:,1:end-to_remove);
%             continue
        else
            I_traces = O.y;
        end
        data = D.y;
        
        [FR, AM, W, Vm, Ri, thr, sr, bump,  pulses, sag_d, a_phase, ap_phase, der_phase, vel_depo, vel_repo, tau_mb, adapt_index, burst_freq, flag] = Analysis_workflow.AP_analysis(experimenter_ID,  data,I_traces, do_plotting, this_exp, p);
%         all_data = Analysis_workflow.AP_analysis(experimenter_ID,  data,I_traces, do_plotting, this_exp, p);
        if flag % this is to avoid adding new data when the recording was not good
           continue
        end
        FR_AP = [FR_AP, FR];
        amp_AP = [amp_AP, AM];
        width_AP = [width_AP, W];
        rest_Vm = [rest_Vm, Vm];
        input_R = [input_R, Ri];
        Firing_threshold = [Firing_threshold, thr];
        SAG_r = [SAG_r, sr]; 
        SAG_d = [SAG_d, sag_d];
        A_bump = [A_bump, bump];
        Vel_depo  = [Vel_depo, vel_depo];
        Vel_repo  = [Vel_repo, vel_repo];
        Tau_mb = [Tau_mb, tau_mb];
        Adapt_index  = [Adapt_index ,adapt_index];
        Burst_freq = [Burst_freq, burst_freq];
        Area_phase = [Area_phase,a_phase];
        AP_phase = [AP_phase,ap_phase];
        DER_phase = [DER_phase,der_phase];
        names(i_exp) = {this_exp};
    catch 
       continue
    end
   
end   

%% Write down to Excel
names_idx = cellfun(@(x) ~isempty(x), names);
NAMES = names(names_idx);
% pulses = GC.inter_pulse_interval.(experimenter_ID) * GC.current_steps.(experimenter_ID);
% current_pulses = GC.inter_pulse_interval.(experimenter_ID) * GC.current_steps.(experimenter_ID);

pulses_str = pulses2str(pulses);
% create table with AP firing rate
this_date_table = table(repmat(recording_date, size(FR_AP,2),1), 'VariableNames', {'Date'});
names_table = cell2table(NAMES, 'VariableNames', {'Cell_ID'});
% delete NaN to adapt to the number of sweeps
to_delete = length(~isnan(FR_AP(:,1))) - length(pulses_str);
FR_AP = FR_AP(1:end-to_delete,:);
FR_AP_table = array2table(FR_AP', 'VariableNames',pulses_str);
amp_AP_table = array2table(amp_AP', 'VariableNames', {'Amplitude'});
width_AP_table = array2table(width_AP', 'VariableNames',{'Max 1st AP width'});
Vm_table = array2table(rest_Vm', 'VariableNames',{'Membrane potential'});
Ri_table = array2table(input_R', 'VariableNames',{'Input Resistance (MOhm)'});
Firing_threshold_table = array2table(Firing_threshold', 'VariableNames',{'AP threshold'});
SAG_ratio_table = array2table(SAG_r', 'VariableNames',{'SAG ratio'});
SAG_diff_table = array2table(SAG_d', 'VariableNames',{'SAG diff'});
Bump_table = array2table(A_bump', 'VariableNames',{'AP Bump'});
Phase_area_table = array2table(Area_phase', 'VariableNames',{'Phase Area'});
Vel_depo_table = array2table(Vel_depo', 'VariableNames',{'Vel Depo'});
Vel_repo_table = array2table(Vel_repo', 'VariableNames',{'Vel Repo'});
Tau_mb_table = array2table(Tau_mb', 'VariableNames',{'Tau mb(ms)'});
Adapt_idx = array2table(Adapt_index', 'VariableNames',{'Adapt idx'});
Burst_freq_table = array2table(Burst_freq', 'VariableNames',{'Burst freq'});
% Phase traces Tables
AP_phase_table = array2table(AP_phase, 'VariableNames', NAMES);
DER_phase_table = array2table(DER_phase, 'VariableNames', NAMES);
DER_phase_table =[repmat(recording_date, length(DER_phase),1), DER_phase_table];
DER_phase_table.Properties.VariableNames(1) = {'date'};
AP_phase_table =[repmat(recording_date, length(AP_phase),1), AP_phase_table];
AP_phase_table.Properties.VariableNames(1) = {'date'};


% sz_FR = size(FR_AP_table,2);
% sz_amp = size(amp_AP_table,1);
% sz_w = size(width_AP_table,1);
%%
filename_xlsx = os.path.join(GC.path_putput_AP_analysis.(experimenter_ID),['AP_frequency_', experimenter_ID, '.xlsx']);


%%
this_table = [this_date_table, names_table,FR_AP_table, amp_AP_table,width_AP_table, Vm_table, Ri_table, Firing_threshold_table, ...
                SAG_ratio_table, SAG_diff_table, Bump_table, Phase_area_table, Vel_depo_table, Vel_repo_table, Tau_mb_table, Adapt_idx, Burst_freq_table];
if exist(filename_xlsx, 'file') && overwrite
    original = readtable(filename_xlsx);
    sz_or = height(original);
    range1 = ['B',char(num2str(sz_or+5))];
    % range for the other sheets
    original2 = readtable(filename_xlsx, 'Sheet', 'Phase AP trace');
    sz_or2 = height(original2);
    range2 =  ['B',char(num2str(sz_or2+5))];
    % Phase DER
    original3 = readtable(filename_xlsx, 'Sheet', 'Phase DER trace');
    sz_or3 = height(original3);
    range3 =  ['B',char(num2str(sz_or3+5))];
       
    % Write table
    writetable(this_table,filename_xlsx,'Sheet',1, 'Range', range1) 
    writetable(AP_phase_table,filename_xlsx,'Sheet','Phase AP trace', 'Range', range2)    
    writetable(DER_phase_table,filename_xlsx,'Sheet','Phase DER trace', 'Range', range3)   
else
    if overwrite
        delete(filename_xlsx)    
        writetable(this_table,filename_xlsx,'Sheet',1, 'Range', 'B1')
        writetable(AP_phase_table,filename_xlsx,'Sheet','Phase AP trace', 'Range', 'B1')
        writetable(DER_phase_table,filename_xlsx,'Sheet','Phase DER trace', 'Range', 'B1')
    else
        writetable(this_table,filename_xlsx,'Sheet',1, 'Range', 'B1')
        writetable(AP_phase_table,filename_xlsx,'Sheet','Phase AP trace', 'Range', 'B1')
        writetable(DER_phase_table,filename_xlsx,'Sheet','Phase DER trace', 'Range', 'B1')

    end
end
DONE = 1;
end

% Helper functions
function puls_str = pulses2str(pulses)
pulses_str = cell(0,0);
for i_p = 1:length(pulses)
    this_pulse = pulses(i_p);
    pulses_str{i_p} = mat2str(this_pulse);
end
puls_str = pulses_str;
end
%%
 %#ok<*AGROW>