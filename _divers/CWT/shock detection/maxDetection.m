function [shockIndexes, thresholdAbsoluteValue] = maxDetection(t, meanWvlt, Freqs, Q, MotherWavelet,...
    ctEdgeEffects, thresholdMode, thresholdValue, maxDetectionMethod, plotMean, meanScale, figName)
%% array size

% array size
if size(t, 1) == 1 && size(meanWvlt, 1) == 1 && size(t, 2) == size(meanWvlt, 2)
    % ok
else
    error('array size problem');
end


%% mean and shock detection

% edge effects
[~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
if length(ctEdgeEffects) == 1
    ctEdgeEffects = [ctEdgeEffects, ctEdgeEffects];
end
TedgeEffects = ctEdgeEffects * DeltaT(Freqs(1));
meanWvlttEdgeEffects = meanWvlt;
meanWvlt((t < t(1)+TedgeEffects(1)) | (t > t(end)-TedgeEffects(2))) = nan;

% mean mean
meanTmeanWvlt = mean(meanWvlt, 'omitnan'); % temporal mean of frequency mean

% shocks
switch thresholdMode
    case 'mean'
        threshold = thresholdValue * meanTmeanWvlt;
    case 'absolute'
        threshold = thresholdValue;
end


shockIndexes = false(1, length(meanWvlt));

switch maxDetectionMethod
    case 'local'
        for k_t = 2:length(meanWvlt)-1
            if isnan(meanWvlt(k_t))
                continue
            end
            if meanWvlt(k_t) >= threshold...
                    && meanWvlt(k_t-1) < meanWvlt(k_t) &&  meanWvlt(k_t) >= meanWvlt(k_t+1)
                shockIndexes(k_t) = true;
            end
        end
        
    case 'global' % max(meanWvlt >= threshold)
        thresholdIndexes = meanWvlt >= threshold;
        
        k1 = 1; k2 = 1;
        while k1 <= length(meanWvlt)
            while k1 <= length(meanWvlt) && ~thresholdIndexes(k1)
                k1 = k1+1;
            end
            k2 = k1;
            while k2 <= length(meanWvlt) && thresholdIndexes(k2)
                k2 = k2+1;
            end
            k2 = k2-1;
            
            [~, km] = max(meanWvlt(k1:k2));
            km = k1-1 + km;
            shockIndexes(km) = true;
            
            k1 = k2+1;
        end
end

% plot
if plotMean
    deltaT = DeltaT(Freqs([1, end]));
    
    plotMeanAndShocks(t, meanWvlt, meanWvlttEdgeEffects, shockIndexes,...
        threshold, meanTmeanWvlt, meanScale, deltaT, figName);
    drawnow;
end

thresholdAbsoluteValue = threshold;

end

