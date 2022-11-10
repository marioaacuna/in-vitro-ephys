function varargout = Mini_EPSC_analysis(experimenter_ID,V_traces, do_plotting, name, p)
% Parse inputs
SR = p{1};
global GC

warning('off')


%% Get the peaks
% maybe it's good to analyse peaks by taking chunks of n datapoints
pts = SR / 100;
length_chunk = round(length(V_traces) / pts);
% check later for an actual length that applies for many types of data
data_reshaped = vec2mat(V_traces, length_chunk);

%% Apply a modified version of https://www.mathworks.com/matlabcentral/fileexchange/19380-postsynaptic-potential-detector
% set parameters (this can be modified later to be upstream)
parameters = struct('minAmp', (-3),... % minimum allowable amplitude for alpha functions (in units of samples)
                'maxAmp', (-5000),... % maximum allowable amplitude
                'minTau', (5),... % minimum allowable tau for alpha functions (in units of samples)
                'maxTau', (2000),... %maximum allowable tau
                'minYOffset', (-100),... % minimum allowable yOffset for alpha functions (in units of mV)
                'maxYOffset', (-30),... % maximum allowable yOffset
                'minDecay', (5),... % minimum allowable decay tau
                'maxDecay', (500),... % maximum allowable decay tau
                'derThresh', (4),... % threshold used to determine if the change of derivative is great enough to separate out two alpha functions
                'closestEPSPs', (5),... % second EPSP is not recognized if it is closer than this value to the previous one (in units of samples)
                'errThresh', (0.08),... % threshold for standard error above which a fit with multiple alphas is attempted
                'dataFilterType', 1,... % 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                'derFilterType', 1,... % 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                'dataFilterLength', 80,... % length of data filter
                'derFilterLength', 50,... % length of derivative filter
                'debugging', (0),... % if set to 1 then debugging figures appear
                'dataStart', 1,... % index of first data point
                'forceDisplay', 0,... % forces a graphical output even if other outputs are taken
                'noFit', 0); % turns off best fitting of start time and decay Tau when 1


if ~do_plotting
     pb = CmdLineProgressBar('Progress ...');
end
%% run iterations per chunk
%  set a counter for events
numEvents = 0;
% create output array for alpha best fit parameters
params = zeros(1,4);
traces = [];
for ichunk =  60: size(data_reshaped,1)
    data = data_reshaped(ichunk,:);
    data_d = detrend(data,1);

    % signal = V_traces;                             % Signal Vector
    % signal = detrend(signal,1);
    L = length(data_d);
    Fs = SR;                                              % Sampling Frequency (Hz)
    % Ts = 1/Fs;                                              % Sampling Interval (sec)
    tv = linspace(0, L/Fs, L);


    [maxValue,indexMax] = max(abs(fft(data_d-mean(data_d))));
    frequency = indexMax * Fs / L;
    data_s = smooth(data_d, 20);
    % filt_data = lowpass(data_s,frequency,Fs);
    filt_data = data_s- median(data_s);
    filt_data = bandstop(filt_data, [frequency/3 1500], Fs);
    filt_data = highpass(filt_data, 4, Fs);
    filt_data = highpass(filt_data, 0.5, Fs);
    
    %% 
    PSPsDown = 1;
    sortedData = sort(filt_data);
    if PSPsDown
        dataMean = mean(sortedData(size(filt_data, 1) * .7:size(filt_data, 1) * .9)); % mean of top 10-30% of data
        dataDev = std(sortedData(size(filt_data, 1) * .6:size(filt_data, 1) * .8)); % standard deviation of top 20-40% of data
        % or noise = data - dataFilt; dataDev = std(noise(1:1000));
        % dataDev = std(highPass(data)) * 3;
    else
        dataMean = mean(sortedData(size(filt_data, 2) * .1:size(filt_data, 2) * .3)); % mean of bottom 10-30% of data
        dataDev = std(sortedData(size(filt_data, 2) * .2:size(filt_data, 2) * .4)); % standard deviation of bottom 20-40% of data
        %dataDev = std(highPass(data)) * 3;
    end
    parameters.errThresh = parameters.errThresh / dataDev;
    % filter the raw data
    parameters.dataFilterType = 0;
    switch parameters.dataFilterType
        case 0 % no filtering
            dataFilt = filt_data;
        case 1 % window filtering
            dataFilt = movingAverage(filt_data, parameters.dataFilterLength);
        case 2 % median filtering
            dataFilt = medfilt1(filt_data, parameters.dataFilterLength);
        case 3 % savitsky-golay
            dataFilt = sgolayfilt(filt_data, 2, parameters.dataFilterLength);
    end

    % filter the derivative
    tempDiff = diff(dataFilt);
    parameters.derFilterType = 2;
    switch parameters.derFilterType
        case 0 % no filtering
            dataDerFilt = tempDiff;
        case 1 % window filtering
            dataDerFilt = movingAverage(tempDiff, parameters.derFilterLength);
        case 2 % median filtering
            dataDerFilt = medfilt1(tempDiff, parameters.derFilterLength);
        case 3 % savitsky-golay
            dataDerFilt = sgolayfilt(tempDiff, 2, parameters.derFilterLength);
    end
    clear tempDiff;
    % filter as per Cohen and Miles 2000
    outData = zeros(size(dataFilt));
   
    if PSPsDown
        for index = 2:length(dataFilt)
            if dataDerFilt(index - 1) < 0
                outData(index) = outData(index - 1) + dataDerFilt(index - 1);
            end
        end
    else
        for index = 2:length(dataFilt)
            if dataDerFilt(index - 1) > 0
                outData(index) = outData(index - 1) + dataDerFilt(index - 1);
            end
        end
    end
    if PSPsDown
        % find where derivative of this function is changing from negative to positive
        functionDer = diff(outData);
        peaks = find((functionDer(2:length(functionDer)) ./ functionDer(1:length(functionDer) -1) < 0 | functionDer(2:length(functionDer)) == 0) & functionDer(1:length(functionDer) - 1) < 0);
    else
        % find where derivative of this function is changing from positive to negative
        functionDer = diff(outData);
        peaks = find((functionDer(2:length(functionDer)) ./ functionDer(1:length(functionDer) -1) < 0 | functionDer(2:length(functionDer)) == 0) & functionDer(1:length(functionDer) - 1) > 0);
    end

    % for each such value greater than derThresh find where the function last
    % began to deviate from 0 and call that an event start
    numStarts = 0;
    whereNull = find(outData == 0);
    whereStarts = zeros(length(peaks), 1);
    wherePeaks = whereStarts;
    for index = 1:length(peaks)
        if abs(outData(peaks(index))) > parameters.derThresh
            numStarts = numStarts + 1;
            whereStarts(numStarts) = whereNull(find(whereNull < peaks(index), 1, 'last'));
            wherePeaks(numStarts) = peaks(index);
        end
    end
    whereStarts(numStarts + 1) = length(outData);
    whereStarts(numStarts + 2:end) = [];
    wherePeaks(numStarts + 1:end) = [];

    %% debug
    if parameters.debugging
        figure
        plot(filt_data);
        line(1:length(outData), outData + mean(dataFilt), 'Color', [1 0 0]);
        %     line(1:length(functionDer), functionDer + mean(dataFilt), 'Color', [0 1 0]);
        line(wherePeaks, dataFilt(wherePeaks), 'Color', [1 0 1], 'linestyle', 'none', 'marker', '+');
        line(whereStarts, dataFilt(whereStarts), 'Color', [0 0 0], 'linestyle', 'none', 'marker', '+');
        if PSPsDown
%             dataMean = mean(dataFilt);
            line([1 length(dataFilt)], -[parameters.derThresh parameters.derThresh] + mean(dataFilt), 'Color', [1 0 0], 'linestyle', ':');
        else
            line([1 length(dataFilt)], [parameters.derThresh parameters.derThresh] + mean(dataFilt), 'Color', [1 0 0], 'linestyle', ':');
        end
        fitLine = line([0 1], [dataMean dataMean], 'Color', [0 0 0]);
        decayLine = line([0 1], [dataMean dataMean], 'Color', [1 0 1]);
        riseStart = line([0 1], [dataMean dataMean], 'linestyle', 'none', 'marker', 'x', 'Color', [0.5430 0.2695 0.0742]);
        fitProps = annotation('textbox',[.8 .1 .1 .1], 'backgroundcolor', [1 1 1]);
        legend('Data', 'Up-Only Function', 'Derivative of U-O Fn', 'Event Peaks', 'Events Starts', 'Threshold for choosing', 'Fit line', 'Decay Fit', 'Rise Line');
    end

    %% start counter
    windows = 1; % so that we include the whole time
    windowSize = length(data);
    num_these_events = 0;
    wherepeaks_to_plot = zeros(1,1);
    for indexEPSP = 1:length(whereStarts) - 1
        if whereStarts(indexEPSP + 1) - whereStarts(indexEPSP) < parameters.closestEPSPs
            continue % skip this iteration if the next fitting spot is really close
        end

        if any(~(whereStarts(indexEPSP)  < windows  - windowSize | whereStarts(indexEPSP)  > windows +  2 * windowSize)) % skip this iteration if the found event is outside of the given windows + a buffer of windowSize on each side in case the location is imprecise
            % must get tau within +/- 50% of the true value or the fitting
            % algorithm does not converge (norm of step goes below boundary
            % condition)
            tauGuess = max(wherePeaks(indexEPSP) - whereStarts(indexEPSP), parameters.minTau);
            ampGuess = dataFilt(wherePeaks(indexEPSP)) - dataFilt(whereStarts(indexEPSP));
            whereStop = whereStarts(indexEPSP + 1) - whereStarts(indexEPSP);
            if whereStop > 4 * tauGuess
                whereStop = 2 * tauGuess;
            end
            [bestParams, residuals] = nlinfit(whereStarts(indexEPSP):whereStarts(indexEPSP) + whereStop, dataFilt(whereStarts(indexEPSP):whereStarts(indexEPSP) + whereStop)', @alphas, [ampGuess, tauGuess, whereStarts(indexEPSP), dataFilt(whereStarts(indexEPSP))], statset('MaxIter', 1000, 'FunValCheck', 'off'));
            %calculate standard error
            bestStdErr = sqrt(residuals' * residuals) / (length(residuals) - 2) / abs(bestParams(1));
            %show the fit line
            if parameters.debugging
                set(fitLine, 'Xdata', whereStarts(indexEPSP):whereStarts(indexEPSP) + whereStop, 'Ydata', alphas(bestParams, whereStarts(indexEPSP):whereStarts(indexEPSP) + whereStop));
                xBounds = get(gca, 'Xlim');
                if xBounds(2) < whereStarts(indexEPSP) + whereStop
                    set(gca, 'Xlim', [whereStarts(indexEPSP) - 50 whereStarts(indexEPSP) + 250]);
                end
            end

            if ~(abs(parameters.minAmp) < abs(bestParams(1))  && abs(parameters.maxAmp * 2.7183) > abs(bestParams(1))  &&... % is the amplitude ok?
                    parameters.minTau < bestParams(2)  && parameters.maxTau > bestParams(2)  &&... % is the tau ok?
                    whereStarts(indexEPSP + 1) - bestParams(3) > bestParams(2) * 0.5 &&... % do we really have enough to fit a function to it?
                    parameters.minYOffset < bestParams(4) && abs(parameters.maxYOffset) > bestParams(4)) ||...% is it in the right place?
                    ~any(~(bestParams(3)  < windows  | bestParams(3)  > windows + windowSize)) % is it in the window that we're looking in
                bestStdErr = inf; % set high so that any other fit is accepted since this one has parameters out of range

            end
            % refit variables
            if bestStdErr < parameters.errThresh

                bestParams(1) = dataFilt(wherePeaks(indexEPSP)) - dataFilt(whereStarts(indexEPSP));
                bestParams(4) = dataFilt(whereStarts(indexEPSP));
                bestParams(2) = wherePeaks(indexEPSP) - bestParams(3); %whereStarts(indexEPSP);
                if ~parameters.noFit
                    % change bestParams(4) to the decay tau
                    whichData = round(bestParams(3) + bestParams(2) * 1.5:min([bestParams(3) + bestParams(2) * 1.5 + 150 whereStarts(indexEPSP + 1)]));
                    if length(whichData) > 5
                        [bestParams(4) FittedDecay] = fitDecaySingle(dataFilt(whichData));
                        residuals = FittedDecay' - dataFilt(whichData);
                        if sqrt(residuals' * residuals) / (length(whichData) - 2) / abs(bestParams(1)) > parameters.errThresh
                            bestParams(4) = NaN;
                            FittedDecay = NaN;
                        end
                    else
                        bestParams(4) = NaN;
                        FittedDecay = NaN;
                    end
                    if isnan(bestParams(4))
                        % make a manual guess
                        if PSPsDown == 0
                            tempDecay = find(dataFilt(wherePeaks(indexEPSP) + 2:min([wherePeaks(indexEPSP) + 90 whereStarts(indexEPSP + 1)])) < (dataFilt(wherePeaks(indexEPSP)) - log(2) * (dataFilt(wherePeaks(indexEPSP)) - min(dataFilt(wherePeaks(indexEPSP) + 2:min([wherePeaks(indexEPSP) + 90 whereStarts(indexEPSP + 1)]))))), 1, 'first');
                        else
                            tempDecay = find(dataFilt(wherePeaks(indexEPSP) + 2:min([wherePeaks(indexEPSP) + 90 whereStarts(indexEPSP + 1)])) > (dataFilt(wherePeaks(indexEPSP)) - log(2) * (dataFilt(wherePeaks(indexEPSP)) - max(dataFilt(wherePeaks(indexEPSP) + 2:min([wherePeaks(indexEPSP) + 90 whereStarts(indexEPSP + 1)]))))), 1, 'first');
                        end
                        if ~isempty(tempDecay)
                            bestParams(4) = tempDecay;
                            FittedDecay = ones([1 length(whichData)]) * dataFilt(wherePeaks(indexEPSP));
                        else
                            bestParams(4) = parameters.minDecay - 1;
                        end
                    end
                    if parameters.debugging
                        set(fitProps, 'String',  [{['\bf Amp: ' num2str(bestParams(1))]};{[' Rise: ' num2str(bestParams(2))]};{[' Start: ' num2str(bestParams(3))]};{[' Decay: ' num2str(bestParams(4))]}]);
                        set(decayLine, 'Xdata', whichData, 'Ydata', FittedDecay);
                        set(gca, 'ylimmode', 'auto');
                        set(gca, 'xlim', [max([1 whereStarts(indexEPSP) - 10]) min([max([1 whereStarts(indexEPSP) - 10]) + max([bestParams(2) * 5 200]) length(data)])]);
                        savedLims = get(gca, 'ylim');
                        set(riseStart, 'xdata', bestParams(3), 'ydata', dataFilt(round(bestParams(3))));
                        set(gca, 'ylim', savedLims);
                        pause
                    end
                else
                    bestParams(4) = parameters.minDecay + 1;
                end

                if bestParams(4) >= parameters.minDecay && bestParams(4) <= parameters.maxDecay &&...
                        abs(bestParams(1)) >= abs(parameters.minAmp) && abs(bestParams(1)) <= abs(parameters.maxAmp)

                    params(numEvents + 1, :) = bestParams;
                    params(numEvents+1,3) = wherePeaks(indexEPSP);
                    wherepeaks_to_plot(num_these_events +1,1) = wherePeaks(indexEPSP);
                    if ~parameters.noFit
                        decayFits{numEvents + 1} = FittedDecay;
                    end
                    
                    % get the event
                     if wherePeaks(indexEPSP) +   100 > windowSize
                        event_trace = NaN(181,1);
                        this_tr = data_s(wherePeaks(indexEPSP) - 80 :  windowSize);
                        event_trace(1:length(this_tr),1) = this_tr; 
                         
                     elseif whereStarts(indexEPSP) < 80
                         event_trace = NaN(181,1);
                         event_trace(whereStarts(indexEPSP) : end,1) = data_s(whereStarts(indexEPSP) : whereStarts(indexEPSP)+ 181-whereStarts(indexEPSP) );
                     else

                         event_trace = data_s(wherePeaks(indexEPSP) - 80 : wherePeaks(indexEPSP) +   100);
                    
                     end
                        
                    traces(numEvents+1,:) = event_trace;
                    % counter
                    numEvents = numEvents + 1;
                    num_these_events = num_these_events +1;
                end

            end
        end
    end
    if parameters.noFit
        decayFits = nan(size(params));
        params(:,4) = nan;
    end

    %
    % % show how well our found trace fits the original data
    % if nargout == 0 || parameters.forceDisplay
    %     if max(data) > 50 || min(data) < -100
    %         parentHandle = scope(data, 1:length(data));
    %     else
    %         parentHandle = scopes(data, 1:length(data));
    %     end
    %     set(0, 'currentfigure', get(parentHandle, 'parent'));
    %     set(gcf, 'currentAxes', parentHandle);
    %
    %     try
    %         for i = 1:size(params, 1)
    %             % plot starts
    %             line(parameters.dataStart + params(i, 3), data(round(params(i, 3))), 'Color', [0 0 1], 'linestyle', 'none', 'marker', '+', 'buttondownfcn', ['set(ancestor(gcbo, ''figure''), ''name'', ''Start time = ' sprintf('%6.1f', params(i,3)) ''')']);
    %             % plot amplitudes
    %             line(parameters.dataStart + [params(i, 3) + params(i, 2) params(i, 3) + params(i, 2)], [data(round(params(i, 3) + params(i, 2))) - params(i, 1) data(round(params(i, 3) + params(i, 2)))], 'Color', [0 1 0], 'buttondownfcn', ['set(ancestor(gcbo, ''figure''), ''name'', [''Rise time = ' sprintf('%4.1f', params(i,2)) ', Amplitude = ' sprintf('%7.2f', params(i,1)) '''])']);
    %             if ~parameters.noFit
    %                 % plot decay tau
    %                 line(parameters.dataStart + ((params(i, 3) + params(i, 2) * 1.5:params(i, 3) + params(i, 2) * 1.5 + numel(decayFits{i}) - 1)), decayFits{i}, 'Color', [1 0 0], 'buttondownfcn', ['set(ancestor(gcbo, ''figure''), ''name'', ''Decay tau = ' sprintf('%6.2f', params(i,4)) ''')']);
    %             end
    %         end
    %     catch
    %         % at least return something
    %
    %     end
    % end

    %
    % figure,
    if do_plotting
        plot(dataFilt)
        hold on
        x =  round(wherepeaks_to_plot);
        plot(x, dataFilt(x), 'r*')
        drawnow
        hold off
    else
       
        pb.print(ichunk,size(data_reshaped,1))
    end  
end



disp(['##################',...
    ' good cell : ', name, ' #############'])
varargout{1}  = params;
varargout{2}  = traces;
% varargout{3}  = width_AP;
% varargout{4}  = Vm;
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
% varargout{19} = flag;



