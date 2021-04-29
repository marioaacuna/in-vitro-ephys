function [FR, amp_AP, width_AP, Vm, AP_thr, SAG_R] = AP_analysis(experimenter_ID,data, do_plotting, name, p)
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
columns= size(data,2);
FR = NaN(length(GC.current_steps.(experimenter_ID)),1);
% origin = 1;
% finish_bsl = 800;
% finish_evk = 7*800+5;
% duration = 0.6;
min_distance = 20;
max_width = 300;
AP_threshold = 20;
% Vm_o = nanmean(data(1: finish_bsl-1));
% if size(data,1) > SR + 1 
% %     keyboard
%     return
% end
% if do_plotting
% %     fig1 = figure();
% end
hold on
for i_data  = 1: columns 
    this_data = data(:,i_data);
    median_evoked = median(this_data(finish_bsl:finish_evk));
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
    if i_data == 1
       % Sag ratio
      sag_p = min(this_data(finish_bsl:finish_evk));
      ss = median(this_data(finish_evk - 100:finish_evk));
    end
end
% if do_plotting
% %     close(fig1)
% end
% determine amplitude of first action potential
ap_idx = intersect(find(FR), find(~isnan(FR)));
trace_to_analyse = data(:,ap_idx(1));
[ap_peaks, ~, ap_w] = findpeaks(trace_to_analyse, 'MinPeakHeight', AP_threshold, 'MinPeakDistance', min_distance, 'MaxPeakWidth', max_width);
% amplitude relative to Vm
Vm = nanmean(trace_to_analyse(1: finish_bsl-1));
amp_AP = (abs(Vm) + ap_peaks(1));
% half width of the first AP
width_AP = 1000*(ap_w(1)/SR);

%% AP threshold
pulses = GC.inter_pulse_interval.(experimenter_ID) * GC.current_steps.(experimenter_ID);
% Set up fittype and options.
ft = fittype('smoothingspline');
% ft = fittype('a/(1+exp(-b*x))');

% ft = fittype( 'a/(1+exp(-b*x))', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NearestInterpolant');
% opts.Display = 'Off';
% opts.StartPoint = [0.957166948242946 0.485375648722841];
% Fit model to data.

[fitresult, ~] = fit(pulses(ap_idx(1)-1:ap_idx(end))', FR(ap_idx(1)-1:ap_idx(end)), ft );
AP_thr = unique(arrayfun(@(y)fzero(@(x)fitresult(x)-1,0),FR(ap_idx(1)-1:ap_idx(end))));
hold on
if do_plotting
    plot(pulses(ap_idx(1)-1:ap_idx(end)), FR(ap_idx(1)-1:ap_idx(end)), 'o')
    hold on
    plot(fitresult)
    plot([AP_thr, AP_thr], [1, 1], 'r*')
    title(name)
    pause
    cla
end

%% SAG ratio
SAG_R = (Vm - sag_p) / (Vm - ss);




