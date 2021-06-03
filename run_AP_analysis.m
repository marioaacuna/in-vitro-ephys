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
    case 'Federica'
        data_path = os.path.join(data_path_root, recording_date{1}, 'traces');
    case 'Liselot'
        data_path = os.path.join(data_path_root, recording_date{1}, [recording_date{1},'.data']);
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
is_Amp = cell2mat(cellfun(@(x) sum(ismember(x,str_exptr)) == length(str_exptr) && ~endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
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
A_bump = [];
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
        
        [FR, AM, W, Vm, Ri, thr, sr, bump,  pulses] = Analysis_workflow.AP_analysis(experimenter_ID,  data,I_traces, do_plotting, this_exp, p);
%         all_data = Analysis_workflow.AP_analysis(experimenter_ID,  data,I_traces, do_plotting, this_exp, p);
        FR_AP = [FR_AP, FR];
        amp_AP = [amp_AP, AM];
        width_AP = [width_AP, W];
        rest_Vm = [rest_Vm, Vm];
        input_R = [input_R, Ri];
        Firing_threshold = [Firing_threshold, thr];
        SAG_r = [SAG_r, sr];
        A_bump = [A_bump, bump];
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
Bump_table = array2table(A_bump', 'VariableNames',{'AP Bump'});
% sz_FR = size(FR_AP_table,2);
% sz_amp = size(amp_AP_table,1);
% sz_w = size(width_AP_table,1);
filename_xlsx = os.path.join(GC.path_putput_AP_analysis.(experimenter_ID),'AP_frequency.xlsx');


%%
this_table = [this_date_table, names_table,FR_AP_table, amp_AP_table,width_AP_table, Vm_table, Ri_table, Firing_threshold_table, SAG_ratio_table, Bump_table];
if exist(filename_xlsx, 'file') && ~overwrite
    original = readtable(filename_xlsx);
    sz_or = height(original);
    range1 = ['B',char(num2str(sz_or+5))];

    writetable(this_table,filename_xlsx,'Sheet',1, 'Range', range1)    
    
else
    if overwrite
        delete(filename_xlsx)    
        writetable(this_table,filename_xlsx,'Sheet',1, 'Range', 'B1')
    else
        writetable(this_table,filename_xlsx,'Sheet',1, 'Range', 'B1')    
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