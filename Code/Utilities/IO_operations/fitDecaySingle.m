function [decayTau1 FittedCurve estimates sse] = fitDecaySingle(yData, PSPtype)
% fits tau to ydata that is slowing in x (only determined by the initial offset)
% [decayTau1 FittedCurve estimates] = fitDecaySingle(yData, PSPtype);
if nargin == 1
   if yData(end) >  yData(1)
       PSPtype = -1;
   else
       PSPtype = 1;
   end
end
if size(yData, 1) < size(yData, 2)
    yData = yData';
end
originalLength = length(yData);
if length(yData) > 100000
    yData(fix(length(yData)/ 500) * 500 + 1:end) = [];
    yData = mean(reshape(yData, 500, []))';
    downSampling = 500;
elseif length(yData) > 10000
    yData(fix(length(yData)/ 50) * 50 + 1:end) = [];                
    yData = mean(reshape(yData, 50, []))';
    downSampling = 50;
else
    downSampling = 1;
end
if length(yData) < 6
    decayTau1 = NaN;
    FittedCurve = NaN;
    estimates = NaN;
    sse = NaN;
    return
end
% generate xData
xdata = (0:length(yData) - 1)';
% Call fminsearch with a random starting point.
start_point = [yData(end) yData(1) - yData(end) length(yData) * -.7];
model = @expfun;
estimates = fminsearch(model, start_point, optimset('MaxFunEvals', 1000, 'Display', 'none'));
% check these for realism
if length(yData) < 200 && (abs(estimates(1) - yData(end)) > abs(yData(1) / 2) || (((estimates(3) >= 0) || estimates(2) * PSPtype == -1)))
   % try fitting with a single exponential
   % [decayTau1 FittedCurve] = fitDecayShort(yData, PSPtype);
   
   % try shaving a few points off of the ends
   if length(yData) > 7
       [decayTau1 FittedCurve estimates sse] = fitDecaySingle(yData(2:end - 2), PSPtype);
       %estimates(1) = estimates(1) + xdata(2) - xdata(5);
       if isnan(decayTau1)
           return
       end
   else
       decayTau1 = NaN;
       FittedCurve = NaN;
       estimates = NaN;
       sse = NaN;
       return
   end
else
    % return exponent with greatest magnitude
    decayTau1 = abs(estimates(3));
end
decayTau1 = decayTau1 * downSampling;
estimates(3) = estimates(3) * downSampling;
if nargout == 0
    figure, plot(xdata .* PSPtype, yData)
    if isnan(decayTau1)
        annotation('textbox',[.25 .5 .5 .1], 'backgroundcolor', [1 1 1], 'String', 'Unable to fit', 'edgecolor', 'none', 'fontsize', 24, 'horizontalalignment', 'center', 'verticalalignment', 'middle');
    else
        line(xdata .* PSPtype, estimates(1) + PSPtype * estimates(2) .* exp(xdata./estimates(3)), 'color', [1 0 0]);
        set(gcf, 'numbertitle', 'off', 'name', ['Tau = ' sprintf('%4.2f', decayTau1)]);    
    end
elseif nargout > 1
    xdata = (0:originalLength - 1)';
    FittedCurve = estimates(1) + PSPtype * estimates(2) .* exp(xdata./estimates(3));
end
if nargout > 3
    if downSampling == 1
        sse = expfun(estimates);
    else
        sse = nan;
    end
end
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for: offest + A * exp(-xdata./tau1) + B * exp(-xdata./tau2) - yData
    function sse = expfun(params)
        ErrorVector = (params(1) + PSPtype * params(2) .* exp(xdata./params(3))) - yData;
        sse = ErrorVector' * ErrorVector;
    end
end
