%% Preamble
% This fuinction will set the parameters for mini analysis fromk igor
% files. In this current version, it reads experiments done for Falkowska
% et al. Future releases will contain implementations for more
% experimenters if needed

function DONE =  run_Mini_SP_EPSC_analysis(experimenter_ID, recording_date, do_plotting, varargin)
clc
DONE = 0;
p = varargin;
overwrite = p{3};
% Set global variables
global GC
% Read general_configs
GC = general_configs();
%%
GC.experimenter_ID = experimenter_ID;
data_path_root =  GC.raw_Mini_data_root_path.(experimenter_ID); % This needs to be added on GC
toolboxes_to_use = {'Igor2Matlab'};
toggle_toolbox(toolboxes_to_use, 'on')
switch experimenter_ID
    case {'Federica', 'Fede_setup_1', 'Fede_setup_1_5HT7','Niels'}
        data_path = os.path.join(data_path_root, recording_date{1}, 'traces');
    case {'Liselot', 'Liselot_setup_1'}
        data_path = os.path.join(data_path_root, '2. opto-5HT7', recording_date{1}, [recording_date{1},'.data']);
    case {'Sri', 'Falk_et_at'}
        data_path = os.path.join(data_path_root, recording_date{1});
    case 'Franziska'
        data_path = os.path.join(data_path_root, recording_date{1});
    otherwise
        keyboard 
end
% for some experiments, it's important to see neurons that do not fire
if strcmp(experimenter_ID, 'Franziska')
    take_astrocyte_shape = 1; % will also take info eventho no peaks are gonna be found
else
    take_astrocyte_shape = 0;
end

% data_path = 'M:\Mario\Fede\traces'; % To be changed later
data_dir = dir(data_path);
% Pick animals
ending = '.ibw';
names ={data_dir.name};
files_in_folder = names(endsWith(names, ending));
files_in_folder = natsort(files_in_folder);

% Isolate the files to analyze
% str_exptr = GC.string_file_selection.(experimenter_ID);
% if any(strcmp(experimenter_ID, {'Niels', 'Sri', 'Liselot_setup_1','Fede_setup_1', 'Falk_et_at', 'Franziska'})) 
%    is_Amp = cell2mat(cellfun(@(x) sum(ismember(x,str_exptr)) == length(str_exptr)+1 && ~endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
% else
%     is_Amp = cell2mat(cellfun(@(x) sum(ismember(x,str_exptr)) == length(str_exptr) && ~endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
% end
% files_to_take = files_in_folder(is_Amp);


% this assumes that all the files are gonna be used otherwise uncomment
% from line 44. No need for an outwave
files_to_take = files_in_folder(:);
% is_outwave = cell2mat(cellfun(@(x) sum(ismember(x,str_exptr)) == length(str_exptr) && endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
% is_outwave = cell2mat(cellfun(@(x) endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
% outwaves_files = files_in_folder(is_outwave); 
% outwaves_files = files_in_folder(:); 
% outwave_identifier = GC.set_string_outwave_selection.(experimenter_ID);
%% Run analysis
% get the parameters to analyse
AMP = [];
Freq = [];
Rise = [];
Decay = [];
TRACES = struct();
names = cell(length(files_to_take),1);
for i_exp = 1:length(files_to_take)
    try
        this_exp = files_to_take{i_exp};
        file_to_read = os.path.join(data_path,this_exp);
        D = IBWread(file_to_read); % read Voltage traces
        % this below might not be needed
%         str_V = strsplit(this_exp,str_exptr);% first characters of the voltage filethat need to match to the Current file
%         start_str_v = [str_V{1}, outwave_identifier] ; 
%         outwave_file_start= char(outwaves_files(startsWith(outwaves_files, start_str_v)));
%         % check for correct identifier
%         if isempty(outwave_file_start)
%             start_str_v = [str_V{1}, str_exptr] ; 
%             outwave_file_start= char(outwaves_files(startsWith(outwaves_files, start_str_v)));       
%         end
%         outwave_file = os.path.join(data_path,outwave_file_start);
%         O =  IBWread(outwave_file); % read Current traces
%         % take only data that has current steps
%         if size(O.y, 2) ~= size(D.y, 2)
% %             disp([this_exp, ' did not correspond to I/O curve'])
%             to_remove = size(O.y, 2) - size(D.y, 2);
%             I_traces = O.y(:,1:end-to_remove);
% %             continue
%         else
%             I_traces = O.y;
%         end
        data = D.y;
        
        [data_out, traces] = Analysis_workflow.Mini_EPSC_analysis(experimenter_ID,  data, do_plotting, this_exp, p);
%         all_data = Analysis_workflow.AP_analysis(experimenter_ID,  data,I_traces, do_plotting, this_exp, p);
        
        Freq = [Freq; numel(data_out(:,1)) / p{2}];
        AMP = [AMP;mean(data_out(:,1))];
        Rise = [Rise;mean(data_out(:,2))/ p{1}]; % this might not be good. Fix later
        Decay = [Decay ;mean(data_out(:,4)) / p{1}];
        this_t = strsplit(this_exp, '.');
        this_t = this_t{1};
        TRACES.(['d',char(recording_date{1})]).(this_t) = traces;
        names(i_exp) = {this_t};

    catch ME
        disp(ME.message)
       continue
    end

end   

%% Write down to Excel
names_idx = cellfun(@(x) ~isempty(x), names);
NAMES = names(names_idx);
% pulses = GC.inter_pulse_interval.(experimenter_ID) * GC.current_steps.(experimenter_ID);
% current_pulses = GC.inter_pulse_interval.(experimenter_ID) * GC.current_steps.(experimenter_ID);

% create the table
this_date_table = table(repmat(recording_date, size(NAMES,1),1), 'VariableNames', {'Date'});
names_table = cell2table(NAMES, 'VariableNames', {'Cell_ID'});
% Freq
FR_table = array2table(Freq, 'VariableNames',{'Frequency'});
% amplitude
amp_table = array2table(AMP, 'VariableNames', {'Amplitude'});
% decay
Decay_table = array2table(Decay, 'VariableNames',{'Decay'});
Rise_table = array2table(Rise, 'VariableNames',{'Rise'});
%%
filename_xlsx = os.path.join(GC.path_putput_AP_analysis.(experimenter_ID),['MINIs_', experimenter_ID, '.xlsx']);


%%
this_table = [this_date_table, names_table,FR_table, amp_table,Decay_table, Rise_table];
if exist(filename_xlsx, 'file') && overwrite
    original = readtable(filename_xlsx);
    sz_or = height(original);
    range1 = ['B',char(num2str(sz_or+5))];
           
    % Write table
    writetable(this_table,filename_xlsx,'Sheet',1, 'Range', range1) 
else
    if overwrite
        delete(filename_xlsx)    
        writetable(this_table,filename_xlsx,'Sheet',1, 'Range', 'B1')
      
    else
        writetable(this_table,filename_xlsx,'Sheet',1, 'Range', 'B1')
    end
end

%% save traces as .mat file
matfilename =  os.path.join(GC.path_putput_AP_analysis.(experimenter_ID),['MINIs_', experimenter_ID, '.mat']);
% check if exists
does_exist = exist(matfilename, "file");
if does_exist
    TRACES_old = load_variable(matfilename, 'TRACES');
    % add this date
    TRACES_old.(['d',char(recording_date{1})]) = TRACES.(['d',char(recording_date{1})]);
    TRACES = TRACES_old;
    
end
save(matfilename, 'TRACES')
    
%% end    
DONE = 1;
end

% % Helper functions
% function puls_str = pulses2str(pulses)
% pulses_str = cell(0,0);
% for i_p = 1:length(pulses)
%     this_pulse = pulses(i_p);
%     pulses_str{i_p} = mat2str(this_pulse);
% end
% puls_str = pulses_str;
% end
%%
 %#ok<*AGROW>