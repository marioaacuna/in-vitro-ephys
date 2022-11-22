function output = alphas(Params, x)
% see if we were passed an array of alphas or just one
if size(Params, 1) > 1
    output = zeros(1, length(x)) + mean(Params(:,4));
    yOffset = 0;
    for index = 1:size(Params, 1)
        Amp = Params(index,1);
        tau = Params(index,2);
        xOffset = Params(index,3);
        
        %eliminate pesky errors
        if tau == 0
            tau = 0.00000000001;
        end        
        tempOutput = Amp * (x - xOffset) / tau .* exp(-(x - xOffset) / tau) + yOffset;
        % make any x parts of the alpha function on the opposite side of the baseline from the alpha equal to the baseline
        if Amp > 0
            tempOutput(tempOutput < yOffset) = yOffset;
        else
            tempOutput(tempOutput > yOffset) = yOffset;
        end
        output = output + tempOutput;
    end
else
    Amp = Params(1);
    tau = Params(2);
    xOffset = Params(3);
    yOffset = Params(4);
        %eliminate pesky errors
        if tau == 0
            tau = 0.00000000001;
        end
        output = Amp * (x - xOffset) / tau .* exp(-(x - xOffset) / tau) + yOffset;
        % make any x parts of the alpha function on the opposite side of the baseline from the alpha equal to the baseline
        if Amp > 0
            output(output < yOffset) = yOffset;
        else
            output(output > yOffset) = yOffset;
        end
end
end
