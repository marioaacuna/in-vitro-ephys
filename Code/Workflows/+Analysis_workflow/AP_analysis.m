function varargout = AP_analysis(experimenter_ID,V_traces, I_traces, do_plotting, name, p, take_astrocyte)
% Parse inputs
SR = p{1};
step_duration = p{4} / 1000;
finish_bsl = SR * (p{3} / 1000);
finish_evk = finish_bsl + (SR *(p{4} /1000));
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
        pulses(i_p) = this_pulse_points(this_pulse_points ~= 0);
    else
        pulses(i_p) = 0;
    end
end
% end
this_pulse_start = 1000* (find(I_traces(:,1) ~= 0,1, 'first') -1) / SR; % start pulse in ms
this_duration_step = round(sum(I_traces(:,i_p) == this_pulse_points(2))) -1;
duration_ms = 1000 * this_duration_step / (size(I_traces,1) / p{2});

if current_steps ~= p{5} || p{6} ~= unique(diff(pulses)) || p{1} ~= (size(I_traces,1) / p{2}) || p{4} ~= duration_ms || p{3} ~= this_pulse_start
    disp (['Experiment: ', name, ' has different number of sweeps or different frequency, PLEASE SELECT DIFFERENT PARAMETERS'])
   
    return
end

%% Calculation Vm
Vm = mean(median(V_traces(1:finish_bsl-1,[1: n_steps])));

%% Setting parameters for peak detection
FR = NaN(length(GC.current_steps.(experimenter_ID)),1);
n_peaks = FR;
% origin = 1;
% finish_bsl = 800;
% finish_evk = 7*800+5;
% duration = 0.6;
if ~take_astrocyte
    min_distance = 20;
    max_width = 300;
    % AP_threshold = 20;
    AP_threshold = 30;
    mpkw = 2;
    min_peak_hight = 5; % for determining single AP afterwards
else
    min_distance = 20;
    max_width = 300;
    % AP_threshold = 20;
    AP_threshold = 7.5;% this is a compromise
    mpkw = 5;
    min_peak_hight = 0;
end

%%
hold on
R_ISI9_ISI1 = NaN; % In case we cannot find more than 10 AP in our traces
do_classify = 1;% plot_adapt = 1;
for i_data  = 1: n_steps 
    this_data = V_traces(:,i_data);
    this_pulse = pulses(i_data);
%     vm_this_trace = median(this_data(1:finish_bsl -10,:));
    median_evoked = median(this_data(finish_bsl:finish_evk));
%     max_evoked = max(this_data(finish_bsl:finish_bsl+20));
%     if max_evoked - vm_this_trace > 0, continue, end
    evoked_trace_norm = this_data(finish_bsl : finish_evk)- median_evoked;
    [peaks, loc, ~] = findpeaks(evoked_trace_norm, 'MinPeakHeight', AP_threshold, 'MinPeakDistance', min_distance, 'MaxPeakWidth', max_width, 'MinPeakWidth', mpkw);
% %     [peaks, loc, w] = findpeaks(this_data,'MinPeakProminence',1, 'MinPeakHeight', mad(this_data(finish_bsl:finish_evk))*1.5, 'MinPeakDistance', min_distance, 'MaxPeakWidth', max_width);\
%     [peaks, loc, w] = findpeaks(this_data,'MinPeakProminence',0.5, 'MinPeakHeight', std(this_data(origin:finish_bsl))*2, 'MinPeakDistance', min_distance, 'MaxPeakWidth', max_width);
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
    FR(i_data) = length(loc)/step_duration;
    n_peaks(i_data) =  length(loc);
    if i_data == 1
       % Sag ratio
      sag_p = min(this_data(finish_bsl:finish_evk));
      ss = median(this_data(finish_evk - ((50 * SR) / 1000):finish_evk -1));
      this_step_Vm = median(this_data(1:finish_bsl-1,1));
      abs_Vm = abs(Vm);
      abs_ss = abs(ss);
      Ri = abs(((abs_ss - abs_Vm)*10^-3) / (pulses(i_data)*10^-9)) /1000; % MOhms
    end
    %% Classify the firing pattern
   
    if length(peaks) >= 10 && do_classify
%         keyboard
        % calculate the diff between peaks
        dist_peaks = (diff(loc) / SR);
        time_bin_activity_n_frames = 2;
        peaks_binned = make_tensor({dist_peaks}, time_bin_activity_n_frames, [], 'median');
%         if plot_adapt
%             figure, heatmap((peaks_binned/ max(peaks_binned)), 'Colormap', jet)% hsv
%             caxis([0 max((peaks_binned/ max(peaks_binned)))])
%             title(name)
%             plot_adapt = 0;
%         end
        
        R_ISI9_ISI1 = peaks_binned(1)/peaks_binned(end);
       
%         if R_ISI9_ISI1 < 0.3
%             figure, heatmap((peaks_binned/ max(peaks_binned)), 'Colormap', jet)% hsv
%             caxis([0 max((peaks_binned/ max(peaks_binned)))])
%             title(name)
%         end
%         R_ISI8_ISI1 = peaks_binned(end-1)/peaks_binned(1);
        do_classify = 0;
%         figure, heatmap(1./loc, 'Colormap', summer)
    end
    %% Do burst frequency AP at 200 pA
    if this_pulse == 200 % So far, take traces at 200 pA, with a threshold of 15ms
        % this calculates the frequency of APs in bursts where the ISI is
        % lower than a given threshold in sec. NaN values are given if there
        % are no spikes at 200pA or if there's no burst.
        ISI_burst_threshold = 0.015;% in sec
        if length(loc) > 1
%             keyboard
            % Determine how much the neuron bursts
            dist_peaks = diff(loc / SR);
            n_peaks_burst = 1 + sum(dist_peaks <= ISI_burst_threshold); % calculate nr of peaks in this time
            burst_freq  = n_peaks_burst/sum(dist_peaks(dist_peaks <= ISI_burst_threshold));
            if isinf(burst_freq), burst_freq = NaN;end %No bursting lower than time threshold 
        else
            burst_freq = NaN;
        end
    end
end
% if do_plotting
% %     close(fig1)
% end
% determine amplitude of first action potential
ap_idx = intersect(find(FR), find(~isnan(FR)));
ap_idx = intersect(ap_idx, find(pulses>0)); % find pulses that are larger than 0 pA to detect an actual AP

if isempty(ap_idx) && ~take_astrocyte
    disp([name, ' is most probably an ASTROCYTE'])
%     amp_AP = NaN;width_AP=NaN; Vm=NaN; AP_thr=NaN; SAG_R=NaN; pulses = NaN;
    flag = 1;
    return
elseif isempty(ap_idx) && take_astrocyte
    disp([name, ' Looks like an ASTROCYTE'])
    flag = 0;
    varargout{1}  = FR;
    varargout{2}  = NaN;
    varargout{3}  = NaN;
    varargout{4}  = Vm;
    varargout{5}  = Ri;
    varargout{6}  = NaN;
    varargout{7}  = NaN;
    varargout{8}  = NaN;
    varargout{9}  = NaN;
    varargout{10} = NaN;
    varargout{11} = NaN;% area
    varargout{12} = NaN(100,1);% phase
    varargout{13} = NaN(100,1);% der
    varargout{14} = NaN;
    varargout{15} = NaN;
    varargout{16} = NaN;
    varargout{17} = R_ISI9_ISI1;
    varargout{18} = burst_freq;
    varargout{19} = flag;

%     amp_AP = NaN;width_AP=NaN; Vm=NaN; AP_thr=NaN; SAG_R=NaN; pulses = NaN;
   
    return
else
    flag = 0;
end


%% Valuaes of 1st AP
trace_to_analyse = V_traces(:,ap_idx(1));
median_evoked_first_AP = median(trace_to_analyse(finish_bsl:finish_evk));
evoked_trace_norm_first_AP = trace_to_analyse(finish_bsl : finish_evk)- median_evoked_first_AP;
[ap_peaks, ap_locs, ap_w] = findpeaks(evoked_trace_norm_first_AP, 'MinPeakHeight', AP_threshold, 'MinPeakDistance', min_distance, 'MaxPeakWidth', max_width);

% amplitude relative to Vm
amp_AP = (ap_peaks(1)+ abs(Vm) + median_evoked_first_AP);

% half width of the first AP
width_AP = 1000*(ap_w(1)/SR);

%% Little shoulder after AP
% smoothing_width = floor(1000 / 1000 * 8 / 1) * 2 + 1;
% Smooth fluorescence trace and compute its first-order derivative
% [~, g] = sgolay(1, smoothing_width);
% Smooth the data
% diff_filter = 1;
% F_smoothed = conv(trace_to_analyse, factorial(diff_filter-1)/(-(1/SR))^(diff_filter-1) * g(:,diff_filter), 'same');
% Apply 1st order derivative to smoothed data
% diff_filter = 2;
% F_derivative = conv(trace_to_analyse, factorial(diff_filter-1)/(-(1/SR))^(diff_filter-1) * g(:,diff_filter), 'same');
AP_loc = ap_locs(1);
evoked_trace_norm_first_AP_smoothed = smooth(evoked_trace_norm_first_AP, SR * 0.001);

if length(ap_locs) > 1
    F_derivative = diff(evoked_trace_norm_first_AP_smoothed(AP_loc+2:AP_loc+ (0.06 * SR)));% take 50 ms after the first AP %% 350 / (1/SR) ; % +250 iniially.
elseif length(ap_locs) == 1
     F_derivative = diff(evoked_trace_norm_first_AP_smoothed(AP_loc+2:end));
end
    
F_der_smoothed= smooth(F_derivative,8);
F_derivative = F_der_smoothed;
% normF_der =(F_derivative / F_derivative);%; * max(trace_to_analyse);
% p1 = 
% p2 = evoked_trace_norm_first_AP(find(F_derivative(AP_loc:end) < 0,100, 'first'))
% pt_of_interest = evoked_trace_norm_first_AP(F_derivative(AP_loc:AP_loc +200) > 0,1);


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
% Keep transitions from falling to rising
% figure, 
% plot(F_derivative)
% for itra = 1:length(zx)
%    hold on
%    plot([zx(itra), zx(itra)], [-0.2 0.2], 'r--')
%    
% end
if length(zx) <3 % in case there are only 2 zero crossings (it could happend if the derivative does not return to 0 within the time window)
    A_bump = NaN; % most probably it will find a second AP instead of a real bump
else
    
    points = [ap_locs(1) + finish_bsl + zx(1) , ap_locs(1) + finish_bsl + zx(2), ap_locs(1) + finish_bsl + zx(3)];
    if length(ap_locs) ~= 1 && length(zx)>2 % In case a sweep has only one action potential, calculate bump directly
        if AP_loc +  zx(3) > ap_locs(2) || (zx(2) - zx(1)) < 5 % if the distance of the third intersection is larger than the position of the second AP, give a NaN; or if the difference between 0-crossings is lower than 5 points
            A_bump = NaN;
        else
            bump = trace_to_analyse(points(1): points(3));
            A_bump = abs(abs(max(bump)) - abs(min(bump)));
            
        end

    else
        if (zx(2) - zx(1)) < 5 % If the distnace of the base (first crossing) and the second crossing is less than 5 points
            A_bump = NaN;
        else
            bump = trace_to_analyse(points(1): points(3));
            A_bump = abs(abs(max(bump)) - abs(min(bump)));
        end
    end
end

%% Determine AP phase plot

% Smooth the long trace
% smoother_trace_to_analyse = smooth(trace_to_analyse, SR * 0.001);
% Find peaks of the long one
[~, ap_locs_long, ~] = findpeaks(trace_to_analyse, 'MinPeakHeight', min_peak_hight, 'MinPeakDistance', min_distance, 'MaxPeakWidth', max_width);
if isempty(ap_locs_long) && ~take_astrocyte % if there's no peak detected, delete max_width
    [~, ap_locs_long, ~] = findpeaks(trace_to_analyse, 'MinPeakHeight', min_peak_hight, 'MinPeakDistance', min_distance);
    disp('peak for phase plot re-checked')
elseif isempty(ap_locs_long) && take_astrocyte % meaning that the amplitude is really small
    disp('%% Very small AP%%')
    [~, ap_locs_long] = max(trace_to_analyse);

end


%
first_AP = trace_to_analyse(ap_locs_long(1)-50:ap_locs_long(1)+50 -1);

% downsample to have the same number of points accross trials
% goal_n_points = 1000;
% actual_n_points = length(first_AP);
% 
% trace_new_phase = resample(first_AP,goal_n_points,actual_n_points);    

% Calculate the derivative 
% trace_new_phase = trace_new_phase(2:end);
der_first = gradient(first_AP)./ (1000/SR); % dV/dt (mV/ms); We use gradient better, cuz it does centered-based diff


% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% time_AP = linspace(1000/SR, length(first_AP)*1000 / SR, length(first_AP));
% der_first_raw =  diff(first_AP, 1000/SR);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure
% figure, plot(first_AP(1:length(first_AP)-1),der_first)
%
% figure, plot(fft(der_first), 'r')
%
% Calculate area phase
AP_to_phase = first_AP;
area_phase = polyarea(AP_to_phase,der_first);
% get the max and min of speed (depol and repol phases, respectively)
speed_depo = max(der_first);
speed_repo = abs(min(der_first));
% trace_phase = [AP_to_phase,der_first ]; % X,Y

% Check figure
% figure, plot(AP_to_phase, der_first)
%% CAlculate tau of membrane depo
trace_to_tau = (V_traces(finish_evk:finish_evk+1000,1));
% trace_to_tau = evoked_trace_norm_first_AP(1:150);
y = trace_to_tau(5:end); % Values of voltage from onse stim to first AP
x = 1000*(1:length(y)) ./ SR;
x = x';
p0= [[ones(size(x)), -exp(-x)]\y; 1];
g = ('a-b*exp(-1/c*x)');
f = fit( x, y, g, 'StartPoint', p0 );
tau_membrane = f.c;
%% AP threshold
% pulses = GC.inter_pulse_interval.(experimenter_ID) * GC.current_steps.(experimenter_ID);
% Set up fittype and options.
ft = fittype('smoothingspline');
% ft = fittype('a/(1+exp(-b*x))');

% ft = fittype( 'a/(1+exp(-b*x))', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NearestInterpolant');
% opts.Display = 'Off';
% opts.StartPoint = [0.957166948242946 0.485375648722841];

% Fit model to data.
[fitresult, ~] = fit(pulses(ap_idx(1)-1:ap_idx(end))', n_peaks(ap_idx(1)-1:ap_idx(end)), ft );
AP_thr = unique(arrayfun(@(y)fzero(@(x)fitresult(x)-1,0),n_peaks(ap_idx(1)-1:ap_idx(end))));
if AP_thr < 0
   AP_thr = pulses(ap_idx(1));
   disp('%%% AP thr was set to first AP current step %%%')
end
hold on
if do_plotting
    plot(pulses(ap_idx(1)-1:ap_idx(end)), n_peaks(ap_idx(1)-1:ap_idx(end)), 'o')
    hold on
    plot(fitresult)
    plot([AP_thr, AP_thr], [1, 1], 'r*')
    title(name)
    pause
    cla
end

%% SAG ratio
SAG_R = (this_step_Vm - sag_p) / (this_step_Vm - ss);
SAG_D = (ss - sag_p);
disp(['##################',...
    ' good cell : ', name, ' #############'])
varargout{1}  = FR;
varargout{2}  = amp_AP;
varargout{3}  = width_AP;
varargout{4}  = Vm;
varargout{5}  = Ri;
varargout{6}  = AP_thr;
varargout{7}  = SAG_R;
varargout{8}  = A_bump;
varargout{9}  = pulses;
varargout{10} = SAG_D;
varargout{11} = area_phase;
varargout{12} = AP_to_phase;
varargout{13} = der_first;
varargout{14} = speed_depo;
varargout{15} = speed_repo;
varargout{16} = tau_membrane;
varargout{17} = R_ISI9_ISI1;
varargout{18} = burst_freq;
varargout{19} = flag;



