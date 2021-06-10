function varargout = EPSP_analysis(~,V_traces, I_traces, do_plotting, name, p)
% Parse inputs
SR = p{1};
start_stim = p{2} / 1000;
testpulse_start = p{3} / 1000;
testpulse_duration = p{5} / 1000;
testpulse_amplitude = p{4}; % in pA
total_duration = p{6};
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
for i_p = 1 : current_steps
    this_pulse_points = unique(I_traces(:,i_p));
    is_0 = sum(this_pulse_points) == 0;
    if ~is_0
        pulses(i_p) = this_pulse_points(this_pulse_points ~= 0);
    else
        pulses(i_p) = 0;
    end
end
current_duration = (sum(I_traces == pulses)-1) / SR;
this_total_duration = (size(I_traces,1)) / SR;
if current_duration ~= testpulse_duration || p{4} ~= unique(pulses) || p{1} ~= ((size(I_traces,1)) / total_duration) || this_total_duration ~= total_duration
    disp (['Experiment: ', name, ' has different number of sweeps or different frequency, PLEASE SELECT DIFFERENT PARAMETERS'])
    return
end

%% Set variables to consider for eact sweep
Vm = NaN(n_sweeps,1);
Ri = Vm;
% EPSP parameters
Amp = Vm;
Rise_time = Vm;
Decay_time = Vm;
Slope = Vm;
Width = Vm;
Onset = Vm; % to be calculated
Sweep_ids = Vm;





%%
% Parameters
min_distance = 0.005;
min_width = 10 / SR;
endpoint = 0.1 * SR; % 100 ms

hold on
for i_data  = 1: n_sweeps 
    this_data = V_traces(:,i_data);
    this_data_smoothed = smooth(this_data, SR * 0.001);
    bsl = this_data_smoothed((start_stim*SR) - (0.1*SR) +1  :start_stim*SR - 50);
    median_bsl= median(bsl);
    this_EPSP_short = this_data_smoothed(start_stim*SR +1 : start_stim*SR + endpoint) - median_bsl;
%     data_start_stim = this_EPSP_short(start_stim*SR :end);
    data_start_stim = this_EPSP_short;
    x = (1:length(data_start_stim)) ./ SR;
    [peaks, loc, w] = findpeaks(data_start_stim, x, 'MinPeakProminence', 0.5,  'MinPeakHeight', 0.2, 'MinPeakDistance', min_distance, 'MinPeakWidth', min_width);
    
    if isempty(peaks) % if there's no peak
        continue
    end
    [~, peak_to_take] = max(peaks);
    this_amp = peaks(peak_to_take);
    % Calculate rise time 
    rise_time = 1000 * risetime(data_start_stim(1:ceil(loc(1)*SR)), SR); % It calculates the rise time to the first peak, if found > 1; in ms
    % Calculate decay time
    if length(rise_time) > 1, rise_time = NaN;end
    y = data_start_stim(loc(peak_to_take)*SR:end); % Values of voltages from the peak to the end 
    x_d = x(ceil(loc(peak_to_take) *SR) :end);
    f = fit(x_d',y,'exp1');
    decay_time = abs(f.b); % in ms ?
   
    % calculate Slope
%     above0 = find(data_start_stim > 0,1, 'first');
    percent_i = this_amp(1) * 0.2;
    percent_f = this_amp(1) * 0.8;
    [~,xi] = min(abs(percent_i-data_start_stim(1: loc(peak_to_take)*SR)));
    [~, xf] =  min(abs(percent_f-data_start_stim(1: loc(peak_to_take)*SR)));
    
    dx = xf - xi;
    
%     dx = x(loc(1)*SR) - x(above0);
    slope = unique(gradient([data_start_stim(xi),data_start_stim(xf)]  ,dx/SR*1000));
    
    % Get the width
    this_width = w(peak_to_take);
    if do_plotting
        cla
        plot(x, data_start_stim, 'k')
        hold on
        plot([loc(peak_to_take), loc(peak_to_take)], [peaks(peak_to_take), peaks(peak_to_take)], 'r*')
        plot(f,x_d,y)
        legend off
        title(['sweep ', num2str(i_data)])
        pause(0.1)
    
    end        
    % calculate Ri
    ss = median(this_data_smoothed(testpulse_start * SR + 0.1*SR : testpulse_start * SR + testpulse_duration * SR));
    this_R = abs(((ss - median_bsl) *10^-3) / (testpulse_amplitude*10^-9)) /1000; % in MOhms
    
    % Calculate onset
    % This calculation takes the 4*std from the baseline value of the trace
    % and subtracts it to the stimulus onset.
    % verify this calculation with other datasets, because the arctifact can be different 
    thr = 4 * std(bsl);
    onset_idx = find(data_start_stim > thr, 1, 'first');
    if isempty (onset_idx) % it could be that there is activity before, increasing the std
        onset_idx = xi;
    end
    this_onset = 1000*((start_stim + onset_idx / SR) - start_stim); % ms

    % Output parameters
    Vm(i_data) = median_bsl;
    Ri (i_data)= this_R;
    % EPSP parameters
    Amp(i_data) = this_amp;
    Rise_time(i_data) = rise_time;
    Decay_time(i_data) = decay_time;
    Slope(i_data) = slope;
    Width(i_data) = this_width;
    Onset(i_data) = this_onset; % to be calculated
    Sweep_ids(i_data) = i_data;
    
end

% Vm = mean(all_Vms);
%% Output arguments
disp(['good cell : ', name])
varargout{1} = Sweep_ids;
varargout{2} = Amp;
varargout{3} = Rise_time;
varargout{4} = Decay_time;
varargout{5} = Slope;
varargout{6} = Width;
varargout{7} = Onset;
varargout{8} = Vm;
varargout{9} = Ri;
