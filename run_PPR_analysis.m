function DONE =  run_PPR_analysis(experimenter_ID, recording_date, do_plotting, varargin)
clc
DONE = 0;
p = varargin;
overwrite = p{8};
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
        unique_datapath = 1;
        
    case 'Liselot'
        data_path = os.path.join(data_path_root, '10. Norbert', recording_date{1}, [recording_date{1},'.data']);
        unique_datapath = 1;
    case 'Kristina'
        initial_path = os.path.join(data_path_root, recording_date{1});
        dir_path = dir(initial_path);
        names_path = {dir_path.name};
        folders_to_analyse = names_path(endsWith(names_path, 'KV'));
        data_path = cell(length(folders_to_analyse),1);
        for i_folder = 1:length(folders_to_analyse)
            data_path{i_folder} = os.path.join(initial_path, folders_to_analyse{i_folder});
        end
        
        if length(data_path) > 1
            unique_datapath = 0;
        else
            data_path = data_path{1};
            unique_datapath = 1;
        end
    case 'Niels'
        data_path_root = ['N:\Niels\Igor\'];
        data_path = os.path.join(data_path_root, recording_date{1},'data');
        unique_datapath = 1;
        
    otherwise
        keyboard
end
if unique_datapath
    DATA = run_analysis(data_path, experimenter_ID, p, do_plotting, recording_date);
else
    DATA = table();
    for i_folder = 1:length(folders_to_analyse)
        data = run_analysis(data_path{i_folder}, experimenter_ID, p, do_plotting, recording_date);
        DATA = [DATA; data];
    end
end

%% Write data
disp('saving data')
the_names = unique(DATA.("Date/Cell"));


for i_name = 1:length(the_names)
    this_name = the_names{i_name};
    data_to_take = find(ismemberCellRows(DATA.("Date/Cell"), {this_name}));
    this_DATA = DATA(data_to_take,:); %#ok<FNDSB>
    filename_xlsx = os.path.join(GC.path_output_EPSP_analysis.(experimenter_ID), [this_name,  '.xlsx']);
    data_fieldnames = fieldnames(this_DATA);
    data_fieldnames(ismember(data_fieldnames, {'Date/Cell', 'trials', 'Properties', 'Row', 'Variables'})) = [];
    % Loop through trials
    trials = this_DATA.trials;
    n_trials = length(trials);
    for i_trial = 1:n_trials
        this_trial = trials{i_trial};
        sheet_name = this_trial;
        idx = ismember(this_DATA.trials,{this_trial});
        sw_id = this_DATA.("Sweep ids"){idx} ;
        amp = this_DATA.Amplitude{idx}; 
        wd = this_DATA.Width{idx};
        rt = this_DATA.("Rise time"){idx};
        dt = this_DATA.("Decay time"){idx};
        on = this_DATA.Onset{idx};
        slp = this_DATA.Slope{idx};
        Ri = this_DATA.("Input Resistance (MOhm)"){idx};
        I_hold = this_DATA.("Ihold"){idx};
        Ra = this_DATA.("Access R"){idx}; 
        PPR = this_DATA.("Paired-pulse ratio"){idx}; 
        % convert back to Table
        this_table = array2table([sw_id, PPR, amp, slp, Ri, Ra, I_hold, on, wd, rt, dt], 'VariableNames', data_fieldnames);
        
        %% Write to Excel
        if exist(filename_xlsx, 'file') && ~overwrite
             writetable(this_table,filename_xlsx,'Sheet',sheet_name, 'Range', 'B1')
%             keyboard
%             original = readtable(filename_xlsx);
%             sz_or = height(original);
%             range1 = ['B',char(num2str(sz_or+5))];
%             writetable(this_table,filename_xlsx,'Sheet',sheet_name, 'Range', range1)
            
        elseif ~exist(filename_xlsx, 'file')

            if overwrite
                delete(filename_xlsx)
                writetable(this_table,filename_xlsx,'Sheet',sheet_name, 'Range', 'B1')
            else
                writetable(this_table,filename_xlsx,'Sheet',sheet_name, 'Range', 'B1')
            end
        end
    end
    
    
    
end
%% DONE
disp('Done!')
DONE = 1;
end

%% Run analysis
    function this_table = run_analysis(data_path, experimenter_ID, p, do_plotting, recording_date)
        global GC
        data_dir = dir(data_path);
        % Pick animals
        ending = '.ibw';
        names ={data_dir.name};
        files_in_folder = names(endsWith(names, ending));
        files_in_folder = natsort(files_in_folder);
        
        % Isolate the files to analyze
        str_exptr = GC.string_file_selection.(experimenter_ID);
        if strcmp(experimenter_ID, 'Liselot')
            is_Amp = cell2mat(cellfun(@(x) sum(ismember(x,str_exptr)) == length(str_exptr) && ~endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
        elseif strcmp(experimenter_ID, 'Niels')
            is_Amp = cell2mat(cellfun(@(x) sum(ismember(x,str_exptr)) == length(str_exptr)+1 && ~endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
        elseif strcmp(experimenter_ID, 'Kristina')
            is_Amp = cell2mat(cellfun(@(x) endsWith(x, 'I1.ibw') && ~endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
            str_exptr = 'I1.ibw';
        end
        files_to_take = files_in_folder(is_Amp);
        % is_outwave = cell2mat(cellfun(@(x) sum(ismember(x,str_exptr)) == length(str_exptr) && endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
        is_outwave = cell2mat(cellfun(@(x) endsWith(x, 'outwave.ibw'), files_in_folder, 'UniformOutput', false));
        outwaves_files = files_in_folder(is_outwave);
        outwave_identifier = GC.set_string_outwave_selection.(experimenter_ID);
        %% Run analysis
        sweep_ids = cell(0,0);
        amp_EPSP =  cell(0,0);
        width_EPSP =  cell(0,0);
        I_hold =  cell(0,0);
        input_R =  cell(0,0);
        access_R = cell(0,0);
%         I_h = cell(0,0);
        slope_EPSP =  cell(0,0);
        onset_EPSP =  cell(0,0);
        risetime_EPSP =  cell(0,0);
        decaytime_EPSP =  cell(0,0);
        ppR = cell(0,0);
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
                V_traces = O.y;
                data = D.y;
                
                [sweep_id, amp, rise, decay, slope, wd, on, i_hold, ri, ra, ppr] = Analysis_workflow.PPR_analysis(experimenter_ID,  data, V_traces, do_plotting, this_exp, p);
                % Output for all experiments on this folder
                
                sweep_ids(i_exp) = {sweep_id};
                amp_EPSP(i_exp) = {amp};
                width_EPSP(i_exp) = {wd} ;
                I_hold(i_exp) = {i_hold};
                input_R(i_exp) = {ri};
                slope_EPSP(i_exp) = {slope};
                onset_EPSP(i_exp) = {on};
                risetime_EPSP(i_exp) = {rise};
                decaytime_EPSP(i_exp) = {decay};
                access_R(i_exp) = {ra}; 
                ppR (i_exp) = {ppr};
                names(i_exp) = {this_exp};
            catch
                continue
            end
            
        end
        
        
        %% Write down to Excel
        names_idx = cellfun(@(x) ~isempty(x), names);
        NAMES = names(names_idx);
        
        
        seps = strsplit(data_path, '\');
        this_folder = (seps{end});
        if strcmp(this_folder, 'data') % For niels' data was something weird.
            this_folder = recording_date{1};
        end
        
        if isempty(NAMES)
%             keyboard
            disp(['File/cell: ', this_folder, ' has no trials'])
            this_table = [];
            return
        end
        % pulses = GC.inter_pulse_interval.(experimenter_ID) * GC.current_steps.(experimenter_ID);
        % current_pulses = GC.inter_pulse_interval.(experimenter_ID) * GC.current_steps.(experimenter_ID);
        
        %         pulses_str = pulses2str(pulses);
        % create table with AP firing rate
        % find sweeps where no peaks were found
        no_peaks_idx = cellfun(@(x) isempty(x), sweep_ids);
        this_date_table = table(repmat({this_folder}, size(sweep_ids(~no_peaks_idx),2),1), 'VariableNames', {'Date/Cell'});
        names_table = cell2table(NAMES, 'VariableNames', {'trials'});
        Sweeps_table = array2table(sweep_ids(~no_peaks_idx)', 'VariableNames',{'Sweep ids'});
        Amp_table = array2table(amp_EPSP(~no_peaks_idx)', 'VariableNames', {'Amplitude'});
        
        
        
        Width_table = array2table(width_EPSP(~no_peaks_idx)', 'VariableNames',{'Width'});
        I_hold_table = array2table(I_hold(~no_peaks_idx)', 'VariableNames',{'Ihold'});
        Ri_table = array2table(input_R(~no_peaks_idx)', 'VariableNames',{'Input Resistance (MOhm)'});
        Slope_table = array2table(slope_EPSP(~no_peaks_idx)', 'VariableNames',{'Slope'});
        Onset_table = array2table(onset_EPSP(~no_peaks_idx)', 'VariableNames',{'Onset'});
        Rise_table = array2table(risetime_EPSP(~no_peaks_idx)', 'VariableNames',{'Rise time'});
        Decay_table = array2table(decaytime_EPSP(~no_peaks_idx)', 'VariableNames',{'Decay time'});
        Access_table = array2table(access_R(~no_peaks_idx)', 'VariableNames',{'Access R'});
        PPR_table = array2table(ppR(~no_peaks_idx)', 'VariableNames',{'Paired-pulse ratio'});
        
        
        %%
       
        
%         filename_xlsx = os.path.join(GC.path_putput_AP_analysis.(experimenter_ID), this_folder,  'AP_frequency.xlsx');
        
        
        %%
        this_table = [this_date_table, names_table,Sweeps_table,PPR_table,  Amp_table, Slope_table,Ri_table, Access_table,I_hold_table, Onset_table, Width_table, Rise_table, Decay_table];
        
        
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