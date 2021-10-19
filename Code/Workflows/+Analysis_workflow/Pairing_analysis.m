function varargout = Pairing_analysis(~,V_traces, I_traces, do_plotting, name, p)
% Parse inputs
SR = p{1};
step_duration = p{4} / 1000;
finish_bsl = SR * (p{3} / 1000);
finish_evk = finish_bsl + (SR *(p{4} /1000))*50; % To have around 100ms of AP
global GC
%% Find peaks
warning('off')
% sampling_freq = 8000;
% start_time = 0.100;
% end_time = 0.7;
n_steps= size(V_traces,2);

% Current traces info
current_steps = size(I_traces,2);
%p{5} = nr of sweeps
%p{6} = inter-sweep interval
pulses = NaN(1, current_steps);
for i_p = 1 : current_steps
    this_pulse_points = unique(I_traces(:,i_p));
    is_0 = sum(this_pulse_points) == 0;
    if ~is_0
        pulses(i_p) = max(this_pulse_points);
    else
        pulses(i_p) = 0;
    end
end
this_duration_step = sum(I_traces(:,i_p) == pulses) / 3; % 3 ^=  number of pulses (to be changed in future releases)
duration_ms = round(1000 * this_duration_step) / (size(I_traces,1) / p{2});
duration_ms = round(duration_ms);%
if n_steps ~= p{5} || p{1} ~= (size(I_traces,1) / p{2}) || p{4} ~= duration_ms % p{6} ~= pulses, not taking it because it doesnt matter
    disp (['Experiment: ', name, ' has different number of sweeps or different frequency, PLEASE SELECT DIFFERENT PARAMETERS'])
    return
end
% 
% %% Calculation Vm
% Vm = mean(median(V_traces(1:finish_bsl-1,[1: n_steps])));
% 
%% Setting parameters for peak detection
% FR = NaN(length(GC.current_steps.(experimenter_ID)),1);
% n_peaks = FR;
amp_AP = NaN(n_steps,3);
width_AP = amp_AP;
after_hyp = NaN(n_steps,1);
Vm = NaN(n_steps,1);
% origin = 1;
% finish_bsl = 800;
% finish_evk = 7*800+5;
% duration = 0.6;
min_distance = 20;
max_width = 300;
AP_threshold = 20;

%%
hold on
% do_classify = 1;% plot_adapt = 1;
for i_data  = 1: n_steps 
    this_data = V_traces(:,i_data);
%     this_pulse = pulses(i_data);
%     vm_this_trace = median(this_data(1:finish_bsl -10,:));
    median_evoked = median(this_data(finish_bsl:finish_evk));
%     max_evoked = max(this_data(finish_bsl:finish_bsl+20));
%     if max_evoked - vm_this_trace > 0, continue, end
    evoked_trace_norm = this_data(finish_bsl : finish_evk)- median_evoked;
    [peaks, loc, ap_w] = findpeaks(evoked_trace_norm, 'MinPeakHeight', AP_threshold, 'MinPeakDistance', min_distance, 'MaxPeakWidth', max_width, 'MinPeakWidth', 2);
    if length(peaks)> 3, keyboard, end
    if do_plotting
        cla
        plot(this_data, 'color', 'k', 'LineWidth', 1.5)
        hold on
        plot(finish_bsl+loc, peaks+median_evoked, 'ro', 'MarkerSize', 10)
        title([name,' ', num2str(i_data)], 'Interpreter', 'none')
        grid on
        %         keyboard
        pause
        cla
    end
    amp_AP(i_data,:) = peaks'; % First AP
    width_AP (i_data,:) = 1000*(ap_w'./SR);
    Vm(i_data) = median(this_data(1:finish_bsl));
    after_hyp (i_data) = min(this_data(0.175*SR:0.25*SR)) - Vm(i_data);
   
end


%% Outputs
disp(['good cell : ', name])

% Vararouts
varargout{1}  = amp_AP;
varargout{2}  = width_AP;
varargout{3}  = after_hyp;
varargout{4}  = Vm;
% varargout{5}  = Ri;
% varargout{6}  = AP_thr;
% varargout{7}  = SAG_R;
% varargout{8}  = A_bump;
% varargout{9}  = pulses;
% varargout{10} = SAG_D;
% varargout{11} = area_phase;
% varargout{12} = AP_to_phase;
% varargout{13} = der_first;
% varargout{14} = speed_depo;
% varargout{15} = speed_repo;
% varargout{16} = tau_membrane;
% varargout{17} = R_ISI9_ISI1;
% varargout{18} = burst_freq;



