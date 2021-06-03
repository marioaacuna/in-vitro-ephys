function varargout = AP_analysis(experimenter_ID,V_traces, I_traces, do_plotting, name, p)
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

if current_steps ~= p{5} || p{6} ~= unique(diff(pulses)) || p{1} ~= (size(I_traces,1) / p{2})
    disp (['Experiment: ', name, ' has different number of sweeps or different frequency, PLEASE SELECT DIFFERENT PARAMETERS'])
    return
end

FR = NaN(length(GC.current_steps.(experimenter_ID)),1);
n_peaks = FR;
% origin = 1;
% finish_bsl = 800;
% finish_evk = 7*800+5;
% duration = 0.6;
min_distance = 20;
max_width = 300;
AP_threshold = 20;
% Vm_o = nanmean(V_traces(1: finish_bsl-1));
% if size(V_traces,1) > SR + 1 
% %     keyboard
%     return
% end
% if do_plotting
% %     fig1 = figure();
% end
hold on
for i_data  = 1: n_steps 
    this_data = V_traces(:,i_data);
%     vm_this_trace = median(this_data(1:finish_bsl -10,:));
    median_evoked = median(this_data(finish_bsl:finish_evk));
%     max_evoked = max(this_data(finish_bsl:finish_bsl+20));
%     if max_evoked - vm_this_trace > 0, continue, end
    evoked_trace_norm = this_data(finish_bsl : finish_evk)- median_evoked;
    [peaks, loc, ~] = findpeaks(evoked_trace_norm, 'MinPeakHeight', AP_threshold, 'MinPeakDistance', min_distance, 'MaxPeakWidth', max_width, 'MinPeakWidth', 2);
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
      Ri = abs((ss*10^-3) / (pulses(i_data)*10^-9)) /1000; % MOhms
    end
end
% if do_plotting
% %     close(fig1)
% end
% determine amplitude of first action potential
ap_idx = intersect(find(FR), find(~isnan(FR)));
ap_idx = intersect(ap_idx, find(pulses>0)); % find pulses that are larger than 0 pA to detect an actual AP

if isempty(ap_idx)
    disp([name, ' is most probably an ASTROCYTE'])
%     amp_AP = NaN;width_AP=NaN; Vm=NaN; AP_thr=NaN; SAG_R=NaN; pulses = NaN;
    return
end
trace_to_analyse = V_traces(:,ap_idx(1));
median_evoked_first_AP = median(trace_to_analyse(finish_bsl:finish_evk));
evoked_trace_norm_first_AP = trace_to_analyse(finish_bsl : finish_evk)- median_evoked_first_AP;
[ap_peaks, ap_locs, ap_w] = findpeaks(evoked_trace_norm_first_AP, 'MinPeakHeight', AP_threshold, 'MinPeakDistance', min_distance, 'MaxPeakWidth', max_width);
% amplitude relative to Vm
% Vm = nanmean(trace_to_analyse(1: finish_bsl-1));
Vm = mean(median(V_traces(1:finish_bsl-1,[1: n_steps])));
amp_AP = (ap_peaks(1)+ abs(Vm) + median_evoked_first_AP);
% amp_AP = (abs(Vm) + ap_peaks(1));
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
F_derivative = diff(evoked_trace_norm_first_AP_smoothed(AP_loc+2:AP_loc+350));% / (1/SR) ; % +250 iniially.
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

points = [ap_locs(1) + finish_bsl + zx(1) , ap_locs(1) + finish_bsl + zx(2), ap_locs(1) + finish_bsl + zx(3)];
if length(ap_locs) ~= 1 % In case a sweep has only one action potential, calculate bump directly
    if AP_loc + zx(3) > ap_locs(2) || (zx(2) - zx(1)) < 5 % if the distance of the third intersection is larger than the position of the second AP, give a NaN; or if the difference between 0-crossings is lower than 5 points
        A_bump = NaN;
    else
        bump = trace_to_analyse(points(1): points(3));
        A_bump = abs(abs(max(bump)) - abs(min(bump)));
        
        %     if startsWith(name, 't76'), keyboard, end
        % %     figure, plot(trace_to_analyse)
        %     figure, plot(trace_to_analyse(points(1): points(3)))
        % %     hold on
        % %     plot([AP_loc + zx(1) + finish_bsl, AP_loc + zx(1) + finish_bsl], [-100 40], 'r--')
        %     title (name)
    end
else
    if (zx(2) - zx(1)) < 5 % If the distnace of the base (first crossing) and the second crossing is less than 5 points
        A_bump = NaN;
    else
        bump = trace_to_analyse(points(1): points(3));
        A_bump = abs(abs(max(bump)) - abs(min(bump)));
    end
end





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
SAG_R = (Vm - sag_p) / (Vm - ss);
disp(['good cell : ', name])
varargout{1} = FR;
varargout{2} = amp_AP;
varargout{3} = width_AP;
varargout{4} = Vm;
varargout{5} = Ri;
varargout{6} = AP_thr;
varargout{7} = SAG_R;
varargout{8} = A_bump;
varargout{9} = pulses;

