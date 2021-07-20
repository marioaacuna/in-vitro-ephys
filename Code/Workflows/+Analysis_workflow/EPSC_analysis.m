function varargout = EPSC_analysis(~,I_traces, V_traces, do_plotting, name, p)
% Parse inputs
SR = p{1};
start_stim = p{2} / 1000;
testpulse_start = p{3} / 1000;
testpulse_duration = p{5} / 1000;
testpulse_amplitude = p{4}; % in pA
total_duration = p{6};
win_to_pk = p{7};
global GC
%% Find peaks
warning('off')
% sampling_freq = 8000;
% start_time = 0.100;
% end_time = 0.7;
n_sweeps= size(I_traces,2);

% Voltage traces info
% current_steps = size(V_traces,2);
voltage_steps = size(V_traces,2); % Voltage steps
%p{5} = nr of sweeps
%p{6} = inter-sweep interval
pulses = NaN(1, voltage_steps);
n_pulses = unique(V_traces);
is_paring_trial = length(n_pulses) == 3; % baseline, AP and hyperpol current

if ~is_paring_trial
    for i_p = 1 : voltage_steps
        this_pulse_points = unique(V_traces(:,i_p));
        is_0 = sum(this_pulse_points) == 0;
        if ~is_0
            pulses(i_p) = this_pulse_points(this_pulse_points ~= 0);
        else
            pulses(i_p) = 0;
        end
    end
    current_duration = (sum(V_traces == pulses)-1) / SR;
    this_total_duration = (size(V_traces,1)) / SR;
    if current_duration ~= testpulse_duration || p{4} ~= unique(pulses) || p{1} ~= ((size(V_traces,1)) / total_duration) || this_total_duration ~= total_duration
        disp (['Experiment: ', name, ' has different number of sweeps or different frequency, PLEASE SELECT DIFFERENT PARAMETERS'])
        return
    end
    endpoint = 0.1 * SR; % 100 ms
    min_distance = 0.005;
    min_width = 10 / SR;
    MinPeakProminence = 0.02;
    MinPeakHeight = 0;% 0.5;
else
    keyboard
    % Set different parameters for AP detection 
    disp(['Experiment: ', name, ' is a pairing trial'])
    [~, max_idx] = max(n_pulses);
    start_stim = find(I_traces, n_pulses(max_idx), 'first');
    start_stim = (start_stim(1) - 1)   /SR;
    endpoint = 0.250 * SR; % 250ms
    min_distance = 0.5;
    min_width = 100 / SR;
    MinPeakProminence = 20;
    MinPeakHeight =  20;
    
    % Pasive membrane params
    [~, min_idx] = min(n_pulses);
    length_test_pulse = length(find(ismember(I_traces, n_pulses(min_idx)))) -1;
    testpulse_duration = length_test_pulse / SR;
    testpulse_start = find( ismember(I_traces, n_pulses(min_idx)), 1, 'first');
    testpulse_start = (testpulse_start - 1) / SR;
end


%% Set variables to consider for eact sweep
Vm = NaN(n_sweeps,1);
Ri = Vm;
Ra = Vm;
Ihold = Vm;
% EPSP parameters
Amp = Vm;
Rise_time = Vm;
Decay_time = Vm;
Slope = Vm;
Width = Vm;
Onset = Vm; % to be calculated
Sweep_ids = Vm;





%%
hold on
for i_data  = 1: n_sweeps 
    this_data = I_traces(:,i_data);
    this_data_smoothed = smooth(this_data, SR * 0.001); % do not take, since it shortenes the peak of the transients
    bsl = this_data((start_stim*SR) - (0.25*SR) +1  :start_stim*SR - 50);
    median_bsl= median(bsl);
    this_EPSC_short = this_data(start_stim*SR +1 : start_stim*SR + endpoint) - median_bsl;
%     data_start_stim = this_EPSP_short(start_stim*SR :end);
    data_start_stim = this_EPSC_short;
    x = (1:length(data_start_stim)) ./ SR;
    
    % Detect max value within the first 10ms or given by the user
    x_10 = (win_to_pk / 1000) * SR;
    
    [min_10, ~] = min(data_start_stim(1:x_10)); % work on it later
    [max_10, loc_max10] = max(data_start_stim(1:x_10)); % work on it later
    ismin = abs(min_10) > abs(max_10);
    if ismin
        data_start_stim = -data_start_stim;
       [max_10, loc_max10] = max(data_start_stim(1:x_10));
%          MinPeakProminence = -MinPeakProminence;
%          MinPeakHeight = MinPeakProminence;
         
    end
%     loc_max10 = loc_max10;
    
    % Detect Peaks
    [peaks, loc, w] = findpeaks(data_start_stim, x,  'MinPeakHeight', MinPeakHeight,  'MinPeakWidth', min_width, 'MinPeakDistance', min_distance);
    
    
    if isempty(peaks) || loc(1) > 0.020 % if there's no peak or if the first peak is found later than 20ms
        if ~ismin
            [~, peak_to_decay] = max(max_10);
        else
            [~, peak_to_decay] = min(max_10);
        end
        loc = loc_max10 / SR;
        w = NaN;
    else

        [~, peak_to_decay] = max(peaks);
    end
    
    this_amp = max_10;
%     this_amp = peaks(peak_to_decay);
    % Calculate rise time 
    try % it could be that there are no EPSPs detected and then it might 
        rise_time = 1000 * risetime(data_start_stim(1:ceil(loc_max10(1))), SR); % It calculates the rise time to the first peak (during the fist 10ms), if found > 1; in ms
    catch
%         keyboard
        rise_time = NaN;
    end
    % Calculate decay time
    if length(rise_time) > 1, rise_time = NaN;end
    y = data_start_stim(loc(peak_to_decay)*SR:end); % Values of current from the peak to the end 
    x_d = x(ceil(loc(peak_to_decay) *SR) :end);
    f = fit(x_d',y,'exp1');
    decay_time = abs(f.b); % in ms ?
   
    % calculate Slope
%     above0 = find(data_start_stim > 0,1, 'first');
    percent_i = this_amp(1) * 0.2;
    percent_f = this_amp(1) * 0.8;
    [~,xi] = min(abs(percent_i-data_start_stim(1: loc(peak_to_decay)*SR)));
    [~, xf] =  min(abs(percent_f-data_start_stim(1: loc(peak_to_decay)*SR)));
    
    dx = xf - xi;
    if dx ~= 0 % happenes when the points are very close together when there's no EPSP
%     dx = x(loc(1)*SR) - x(above0);
        slope = unique(gradient([data_start_stim(xi),data_start_stim(xf)]  ,dx/SR*1000));
    else
        slope = NaN;
    end
    % Get the width
    this_width = w(peak_to_decay) * 1000;
    if do_plotting
        cla
        plot(x, -data_start_stim, 'k')
        hold on
        plot(loc_max10/SR, -this_amp,  'r*')
        plot(f,x_d,y)
        legend off
        title([name, ' ; sweep ', num2str(i_data)])
        pause(0.1)
    
    end        
    % calculate Ri % Ra "tau * dV/dQt"
    ss = median(this_data_smoothed(testpulse_start * SR + 0.01*SR : testpulse_start * SR + testpulse_duration * SR -50));
    this_R = ((testpulse_amplitude*10^-3)) / abs(((ss - median_bsl) *10^-9)) / 1000000; % in MOhms
    test_pulse_trace = (this_data(testpulse_start *SR -5  : testpulse_start*SR -1 + 0.05*SR - 1)) - median_bsl; % bring the testpulse to 0
   
    [~, testpulse_loc] = min(test_pulse_trace);
    y_t =  test_pulse_trace(testpulse_loc + 5 :end);
    x_t_ms = (1:length(test_pulse_trace)) ./ SR;
    x_t = x_t_ms(ceil(testpulse_loc + 5 :end));
    f_t = fit(x_t',y_t,'exp2');
    tau_test = abs(f_t.d);
    tau_double = ((f_t.a * f_t.c) + (f_t.c * f_t.d)) / (f_t.a + f_t.c);
    
    F_derivative = gradient(test_pulse_trace);
    % Find 0-crossings in 1st derivative (i.e., sign of product of consecutive
    % samples is negative)
    zx = find(sign(F_derivative(1:end-1).*F_derivative(2:end)) < 0);
    % Remove spurious points
    zx(zx<1) = [];
    zx(zx>=length(F_derivative)) = [];
    % Get the sign of points around 0-crossings
    yx = [F_derivative(zx) F_derivative(zx+1)];
    % Keep transitions from rising to falling
    pos_zx = zx(yx(:,1)>=0 & yx(:,2)<0);
    
    % Find the 2nd 0-crossing
%     Qt = trapz(test_pulse_trace(1: pos_zx(2)));  %f_t.c .* f_t.dpos_zx(2) 
    % Cut the trace to only have the area under the transient
    new_trace = test_pulse_trace(test_pulse_trace <= test_pulse_trace(pos_zx(2))) -(test_pulse_trace(pos_zx(2)));
    Qt = trapz(new_trace);
    this_Ra = abs(tau_double) * abs((testpulse_amplitude ) /Qt);
%     this_Ra = NaN;
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
    
    % Rm = (ΔV - Ra * ΔI) / ΔI
    Ihold(i_data) = median_bsl * 1000; % in pA
    Ri (i_data) = this_R;
    Ra (i_data) = this_Ra;
    % EPSP parameters
    Amp(i_data) = this_amp * 1000; % in pA
    Rise_time(i_data) = rise_time; % in ms
    Decay_time(i_data) = decay_time;  % in ms
    Slope(i_data) = slope; 
    Width(i_data) = this_width; % in ms
    Onset(i_data) = this_onset; % in ms
    Sweep_ids(i_data) = i_data;
    
end



% Vm = mean(all_Vms);
%% Output arguments
disp(['good cell : ', name])
varargout{1}  = Sweep_ids;
varargout{2}  = Amp;
varargout{3}  = Rise_time;
varargout{4}  = Decay_time;
varargout{5}  = Slope;
varargout{6}  = Width;
varargout{7}  = Onset;
varargout{8}  = Ihold;
varargout{9}  = Ri;
varargout{10} = Ra;

