function varargout = EPSP_summation_analysis(~,V_traces, I_traces, do_plotting, name, p)
% Parse inputs
SR = p{1};
start_stim = p{2} / 1000;
testpulse_start = p{3} / 1000;
if ~isempty(p{5})
    testpulse_duration = p{5} / 1000;
end
testpulse_amplitude = p{4}; % in pA
total_duration = p{6};
win_to_pk = p{7};
global GC
%% Find peaks
warning('off')
% sampling_freq = 8000;
% start_time = 0.100;
% end_time = 0.7;
n_sweeps= size(V_traces,2);

% Current traces info
current_steps = size(I_traces,2);
%p{5} = nr of sweeps
%p{6} = inter-sweep interval
pulses = NaN(1, current_steps);
% is_paring_trial = length(n_pulses) == 3; % baseline, AP and hyperpol current


if   p{1} ~= ((size(I_traces,1)) / total_duration) 
    disp (['Experiment: ', name, ' wrong parameters'])
    return
end





% 
% 
% 
% if ~is_paring_trial
%     for i_p = 1 : current_steps
%         this_pulse_points = unique(I_traces(:,i_p));
%         is_0 = sum(this_pulse_points) == 0;
%         if ~is_0
%             pulses(i_p) = this_pulse_points(this_pulse_points ~= 0);
%         else
%             pulses(i_p) = 0;
%         end
%     end
%     current_duration = (sum(I_traces == pulses)-1) / SR;
%     this_total_duration = (size(I_traces,1)) / SR;
%     if current_duration ~= testpulse_duration || p{4} ~= unique(pulses) || p{1} ~= ((size(I_traces,1)) / total_duration) || this_total_duration ~= total_duration
%         disp (['Experiment: ', name, ' has different number of sweeps or different frequency, PLEASE SELECT DIFFERENT PARAMETERS'])
%         return
%     end
%     endpoint = 0.1 * SR; % 100 ms
%     min_distance = 0.005;
%     min_width = 10 / SR;
%     MinPeakProminence = 0.2;
%     MinPeakHeight =  0.5;
% else
%     % Set different parameters for AP detection 
%     disp(['Experiment: ', name, ' is a pairing trial'])
%     [~, max_idx] = max(n_pulses);
%     start_stim = find(I_traces, n_pulses(max_idx), 'first');
%     start_stim = (start_stim(1) - 1)   /SR;
%     endpoint = 0.250 * SR; % 250ms
%     min_distance = 0.005;
%     min_width = 10 / SR;
%     MinPeakProminence = 20;
%     MinPeakHeight =  20;
%     
%     % Pasive membrane params
%     [~, min_idx] = min(n_pulses);
%     length_test_pulse = length(find(ismember(I_traces, n_pulses(min_idx)))) -1;
%     testpulse_duration = length_test_pulse / SR;
%     testpulse_start = find( ismember(I_traces, n_pulses(min_idx)), 1, 'first');
%     testpulse_start = (testpulse_start - 1) / SR;
% end


%% Set variables to consider for eact sweep
Vm = NaN(n_sweeps,1);
Area = Vm;
% EPSP parameters
Amp = Vm;
Ratio = Vm;
Rise_time = Vm;
Decay_time = Vm;
Slope = Vm;
Width = Vm;
Onset = Vm; % to be calculated
Sweep_ids = Vm;

% Peak detection parameters
endpoint = 0.1 * SR; % 100 ms
min_distance = 0.020;
min_width = 0.020;
MinPeakProminence = 1.5;
MinPeakHeight =  1.5;

% detect if this data based on only one epsp or more
this_data = I_traces(:,2);
this_data_smoothed = smooth(this_data, SR * 0.001);
bsl = this_data_smoothed(1 :start_stim*SR-1 );
median_bsl= median(bsl);
this_EPSP_short = this_data_smoothed(start_stim*SR +1 : end) - median_bsl;
%     data_start_stim = this_EPSP_short(start_stim*SR :end);
data_start_stim = this_EPSP_short;
x = (1:length(data_start_stim)) ./ SR;
[I_peaks, I_loc, I_w] = findpeaks(data_start_stim, x, 'MinPeakProminence', 0,...
    'MinPeakHeight',0, 'MinPeakDistance', 0, 'MinPeakWidth', 0);
n_peaks = length(I_peaks);


%% detect the peaks of each stimulation
I_total_peaks = NaN(n_sweeps, n_peaks);
for is = 1:n_sweeps
    
    this_data = I_traces(:,is);
    this_data_smoothed = smooth(this_data, SR * 0.001);
    bsl = this_data_smoothed(1 :start_stim*SR-1 );
    median_bsl= median(bsl);
    this_EPSP_short = this_data_smoothed(start_stim*SR +1 : end) - median_bsl;
    %     data_start_stim = this_EPSP_short(start_stim*SR :end);
    data_start_stim = this_EPSP_short;
    x = (1:length(data_start_stim)) ./ SR;
    [I_peaks, I_loc, I_w] = findpeaks(data_start_stim, x, 'MinPeakProminence', 0,...
        'MinPeakHeight',0, 'MinPeakDistance', 0, 'MinPeakWidth', 0);
    if ~isempty(I_peaks)
        I_total_peaks(is, :) = round(I_peaks)';
    else
         I_total_peaks(is, :) = zeros(1,n_peaks);
    end
    %     n_peaks = length(I_peaks);
     
end

%%
Amp = NaN(n_sweeps, n_peaks);
Ratio =  NaN(n_sweeps,1);
if n_peaks == 1
     is_only_1_protocol = 1;
     Amp = NaN(n_sweeps, n_peaks);
     Ratio = Vm;

elseif n_peaks == 0
%       keyboard
    
else
    
%     keyboard
    is_only_1_protocol = 0;
    n_peaks_this_protocol = length(I_loc);
    peaks_loc = I_loc;
end

%%
if is_only_1_protocol
    
    hold on
    for i_data  = 1: n_sweeps
        
        % Detect Peaks voltage trace
        this_data = V_traces(:,i_data);
        this_data_smoothed = smooth(this_data, SR * 0.001);
        bsl = this_data_smoothed(1 :start_stim*SR-1 );
        median_bsl= median(bsl);
        this_EPSP_short = this_data_smoothed(start_stim*SR +1 : start_stim*SR +1 + 0.1*SR) - median_bsl;
        %     data_start_stim = this_EPSP_short(start_stim*SR :end);
        data_start_stim = this_EPSP_short;
        x = (1:length(data_start_stim)) ./ SR;
        [peaks, loc, w] = findpeaks(data_start_stim, x, 'MinPeakProminence', MinPeakProminence,...
            'MinPeakHeight',MinPeakHeight, 'MinPeakDistance', min_distance, 'MinPeakWidth', min_width);
        
        if isempty(peaks) %|| loc(1) > 0.015 % if there's no peak or if the first peak is found later than 15ms
            this_amp = NaN;
            loc = NaN;
            w = NaN;
            this_area = NaN;
            this_ratio = NaN;
            Sweep_ids(i_data) = i_data;
            continue
            %         [~, peak_to_decay] = max(max_10);
            %         loc = loc_max10 / SR;
            %         w = NaN;
            %     else
            %
            %         [~, peak_to_decay] = max(peaks);
        end
        
        this_amp = peaks;
        
        
        F_derivative = gradient(smooth(this_EPSP_short, 100));
        
        % Get 0-crossings
        zx = find(sign(F_derivative(1:end-1).*F_derivative(2:end)) < 0);
        % Remove spurious points
        zx(zx<1) = [];
        zx(zx>=length(F_derivative)) = [];
        yx = [F_derivative(zx) F_derivative(zx+1)];
        % Keep transitions from rising to falling
        pos_zx = zx(yx(:,1)>=0 & yx(:,2)<0);
        
        if size(pos_zx,1) == 2
            epsp_to_area = this_EPSP_short(1:pos_zx(2)); 
        elseif size(pos_zx,1) >= 3
            epsp_to_area = this_EPSP_short(1:pos_zx(3));
        elseif size(pos_zx,1) == 1
            epsp_to_area = this_EPSP_short(1:end);
        end
        this_area = sum(epsp_to_area); %trapz(epsp_to_area)
        
        if do_plotting
            cla
            plot(x, data_start_stim, 'k')
            hold on
            plot(loc, this_amp,  'r*')
%             plot(f,x_d,y)
            legend off
            title([name, ' ; sweep ', num2str(i_data)])
            pause(0.1)
            
        end
        
        % Output parameters
        Vm(i_data) = median_bsl;
        Area (i_data)= this_area;
        % EPSP parameters
        Amp(i_data) = this_amp;
        Ratio(i_data) = this_amp/this_amp;  
%         Rise_time(i_data) = [];%rise_time;
%         Decay_time(i_data) = [];%decay_time;
%         Slope(i_data) =[]; %slope;
        Width(i_data) = w;
%         Onset(i_data) = [];%this_onset; % to be calculated
        Sweep_ids(i_data) = i_data;
        
    end
    
    disp('only one injection per sweep')
else
    hold on
    for i_data  = 1: n_sweeps
        
        
        % Detect max value within the first 10ms or given by the user
        %     x_10 = (win_to_pk / 1000) * SR;
        %     [max_10, loc_max10] = max(data_start_stim(1:x_10));
        %     loc_max10 = loc_max10;
        
        
        % Detect Peaks voltage trace
        this_data = V_traces(:,i_data);
        this_data_smoothed = smooth(this_data, SR * 0.001);
        bsl = this_data_smoothed(1 :start_stim*SR-1 );
        median_bsl= median(bsl);
        this_EPSP_short = this_data_smoothed(start_stim*SR +1 : end) - median_bsl;
        %     data_start_stim = this_EPSP_short(start_stim*SR :end);
        data_start_stim = this_EPSP_short;
        x = (1:length(data_start_stim)) ./ SR;
        
        % set if there are peaks
        if sum(I_total_peaks(i_data,:)) == 0
            p = [];
        else
            p = 1;
        end
        if ~isempty(p)
            
            x_peaks = I_loc * SR;
            peaks = [];
            loc = [];
            w=[];
            for ip = 1:n_peaks_this_protocol
                [this_max, this_loc] =max(data_start_stim(x_peaks(ip)-(0.0015*SR): x_peaks(ip) + (0.0150 *SR)));
                this_loc = this_loc + x_peaks(ip) - (0.0016*SR);
                peaks(ip) = this_max;
                loc(ip) = this_loc;   
            end
            
        else %if isempty(peaks) %|| loc(1) > 0.015 % if there's no peak or if the first peak is found later than 15ms
            this_amp = NaN(size(I_total_peaks));
            loc = NaN(size(I_total_peaks));
            w = NaN(size(I_total_peaks));
            this_area = NaN;
            Sweep_ids (i_data) = i_data;
            
            
            if do_plotting
                cla
                plot(x, data_start_stim, 'k')
                hold on
                plot(loc, this_amp,  'r*')
%                 plot(f,x_d,y)
                legend off
                title([name, ' ; sweep ', num2str(i_data)])
                pause(0.1)
                
            end
            
            continue
        end        
        % stop if an AP was detected
        is_AP = (peaks > 30);
        if any(is_AP)
            varargout{1} = Sweep_ids;
            varargout{2} = Amp;
            varargout{3} = I_total_peaks;
            varargout{4} = [];
            varargout{5} = [];
            varargout{6} = [];
            varargout{7} = Ratio;
            varargout{8} = Vm;
            varargout{9} = Area;
             disp([num2str(n_peaks_this_protocol),' injections per sweep'])
            return
            continue
        end
        
        this_amp = peaks;
        
        
        F_derivative = gradient(smooth(this_EPSP_short, 100));
        
        % Get 0-crossings
        zx = find(sign(F_derivative(1:end-1).*F_derivative(2:end)) < 0);
        % Remove spurious points
        zx(zx<1) = [];
        zx(zx>=length(F_derivative)) = [];
        yx = [F_derivative(zx) F_derivative(zx+1)];
        % Keep transitions from rising to falling
        pos_zx = zx(yx(:,1)>=0 & yx(:,2)<0);
        idx_to_take = 2*n_peaks_this_protocol -1;
        if idx_to_take > length(pos_zx)
            epsp_to_area = this_EPSP_short(1:pos_zx(end));
        else
            epsp_to_area = this_EPSP_short(1:pos_zx(idx_to_take));
        end
        this_area = sum(epsp_to_area); %trapz(epsp_to_area)
        
        
        
        
        %     this_amp = peaks;
        % %     this_amp = peaks(peak_to_decay);
        %     % Calculate rise time
        %     try % it could be that there are no EPSPs detected and then it might
        %         rise_time = 1000 * risetime(data_start_stim(1:ceil(loc_max10(1))), SR); % It calculates the rise time to the first peak (during the fist 10ms), if found > 1; in ms
        %     catch
        % %         keyboard
        %         rise_time = NaN;
        %     end
        %     % Calculate decay time
        %     if length(rise_time) > 1, rise_time = NaN;end
        %     y = data_start_stim(loc(peak_to_decay)*SR:end); % Values of voltages from the peak to the end
        %     x_d = x(ceil(loc(peak_to_decay) *SR) :end);
        %     f = fit(x_d',y,'exp1');
        %     decay_time = abs(f.b); % in ms ?
        %
        %     % calculate Slope
        % %     above0 = find(data_start_stim > 0,1, 'first');
        %     percent_i = this_amp(1) * 0.2;
        %     percent_f = this_amp(1) * 0.8;
        %     [~,xi] = min(abs(percent_i-data_start_stim(1: loc(peak_to_decay)*SR)));
        %     [~, xf] =  min(abs(percent_f-data_start_stim(1: loc(peak_to_decay)*SR)));
        %
        %     dx = xf - xi;
        %     if dx ~= 0 % happenes when the points are very close together when there's no EPSP
        % %     dx = x(loc(1)*SR) - x(above0);
        %         slope = unique(gradient([data_start_stim(xi),data_start_stim(xf)]  ,dx/SR*1000));
        %     else
        %         slope = NaN;
        %     end
        %     % Get the width
        %     this_width = w(peak_to_decay) * 1000;
        if do_plotting
            cla
            plot(x, data_start_stim, 'k')
            hold on
            plot(loc/SR, this_amp,  'r*')
%             plot(f,x_d,y)
            legend off
            title([name, ' ; sweep ', num2str(i_data)])
            pause(0.1)
            
        end
        %     % calculate Ri
        %     ss = median(this_data_smoothed(testpulse_start * SR + 0.1*SR : testpulse_start * SR + testpulse_duration * SR));
        %     this_R = abs(((ss - median_bsl) *10^-3) / (testpulse_amplitude*10^-9)) /1000; % in MOhms
        %
        %     % Calculate onset
        %     % This calculation takes the 4*std from the baseline value of the trace
        %     % and subtracts it to the stimulus onset.
        %     % verify this calculation with other datasets, because the arctifact can be different
        %     thr = 4 * std(bsl);
        %     onset_idx = find(data_start_stim > thr, 1, 'first');
        %     if isempty (onset_idx) % it could be that there is activity before, increasing the std
        %         onset_idx = xi;
        %     end
        %     this_onset = 1000*((start_stim + onset_idx / SR) - start_stim); % ms
        
        % Output parameters
        Vm(i_data) = median_bsl;
        Area (i_data)= this_area;
        % EPSP parameters
        Amp(i_data,:) = this_amp;
        Ratio(i_data) = this_amp(end) / this_amp(1); 
        Rise_time(i_data) = [];%rise_time;
%         Decay_time(i_data) = [];%decay_time;
%         Slope(i_data) =[]; %slope;
%         Width(i_data) = [];
        Onset(i_data) = [];%this_onset; % to be calculated
        Sweep_ids(i_data) = i_data;
    end
   
end





% Vm = mean(all_Vms);
%% Output arguments
disp(['good cell : ', name])
varargout{1} = Sweep_ids;
varargout{2} = Amp;
varargout{3} = I_total_peaks;
varargout{4} = [];
varargout{5} = [];
varargout{6} = [];
varargout{7} = Ratio;
varargout{8} = Vm;
varargout{9} = Area;

