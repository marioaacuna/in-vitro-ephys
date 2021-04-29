function run_AP_analysis(experimenter_ID, recording_date, do_plotting, varargin)

p = varargin;
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

%% Run analysis
FR_AP = [];%NaN(length(GC.current_steps.(experimenter_ID)),1);
amp_AP = [];
width_AP = [];
rest_Vm = [];
Firing_threshold = [];
SAG_r = [];
% do_plotting = 0;
names = cell(length(files_to_take),1);
for i_exp = 1:length(files_to_take)
    try
        this_exp = files_to_take{i_exp};
        file_to_read = os.path.join(data_path,this_exp);
        D = IBWread(file_to_read);
        data = D.y;
        [FR, AM, W, Vm, thr, sr] = Analysis_workflow.AP_analysis(experimenter_ID,  data, do_plotting, this_exp, p);
        FR_AP = [FR_AP, FR];
        amp_AP = [amp_AP, AM];
        width_AP = [width_AP, W];
        rest_Vm = [rest_Vm, Vm];
        Firing_threshold = [Firing_threshold, thr];
        SAG_r = [SAG_r, sr];
    catch
        continue
    end
    names(i_exp) = {this_exp};
end   

%% Write down to Excel
names_idx = cellfun(@(x) ~isempty(x), names);
NAMES = names(names_idx);
pulses = GC.inter_pulse_interval.(experimenter_ID) * GC.current_steps.(experimenter_ID);
pulses_str = pulses2str(pulses);
% create table with AP firing rate
this_date_table = table(repmat(recording_date, size(FR_AP,2),1), 'VariableNames', {'Date'});
names_table = cell2table(NAMES, 'VariableNames', {'Cell_ID'});
FR_AP_table = array2table(FR_AP', 'VariableNames',pulses_str);
amp_AP_table = array2table(amp_AP', 'VariableNames', {'Amplitude'});
width_AP_table = array2table(width_AP', 'VariableNames',{'Max 1st AP width'});
Vm_table = array2table(rest_Vm', 'VariableNames',{'Membrane potential'});
Firing_threshold_table = array2table(Firing_threshold', 'VariableNames',{'AP threshold'});
SAG_ratio_table = array2table(SAG_r', 'VariableNames',{'SAG ratio'});
% sz_FR = size(FR_AP_table,2);
% sz_amp = size(amp_AP_table,1);
% sz_w = size(width_AP_table,1);
filename_xlsx = os.path.join(GC.path_putput_AP_analysis.(experimenter_ID),'AP_frequency.xlsx');


%%
this_table = [this_date_table, names_table,FR_AP_table, amp_AP_table,width_AP_table, Vm_table, Firing_threshold_table, SAG_ratio_table];
if exist(filename_xlsx, 'file')
    original = readtable(filename_xlsx);
    sz_or = height(original);
    range1 = ['B',char(num2str(sz_or+2))];

    writetable(this_table,filename_xlsx,'Sheet',1, 'Range', range1)    
    
else
    writetable(this_table,filename_xlsx,'Sheet',1, 'Range', 'B1')    
end
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