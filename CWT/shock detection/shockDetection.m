function spectrums = shockDetection(t, X, freqsMean, freqsSpectrum, QMean, QSpectrum, MotherWaveletMean,...
    MotherWaveletSpectrum, ctEdgeEffectsMean, ZeroPadding,...
    meanFunc, thresholdMode, thresholdValue, maxDetectionMethod, varargin)
%CHOCDETECTION Summary of this function goes here
%   Wavelet : CWT(i,j)  i freq index & j time index

p = inputParser;
% parametres par defaut

% meanFunctionDef = @(x) abs(x).^2;
% thresholdModeDef = 'mean'; % 'mean', 'absolute'
% maxDetectionMethodDef = 'local'; % 'local', 'global'  (detect local maxima or global maxima)
plotMeanDef = true;
plotSpectrumDef = true;
meanScaleDef = 'lin';
spectrumFrequencyScaleDef = nan;
spectrumScaleDef = 'lin';

addParameter(p, 'plotMean', plotMeanDef);
addParameter(p, 'plotSpectrum', plotSpectrumDef);
addParameter(p, 'meanScale', meanScaleDef);
addParameter(p, 'spectrumFrequencyScale', spectrumFrequencyScaleDef);
addParameter(p, 'spectrumScale', spectrumScaleDef);

parse(p, varargin{:});

plotMean = p.Results.plotMean;
plotSpectrum = p.Results.plotSpectrum;
meanScale = p.Results.meanScale;
spectrumFrequencyScale = p.Results.spectrumFrequencyScale;
spectrumScale = p.Results.spectrumScale;

if any(isnan(spectrumFrequencyScale))
    if max(abs(diff(freqsSpectrum) / mean(diff(freqsSpectrum)) - 1)) < 1e-3
        spectrumFrequencyScale = 'lin';
    elseif max(abs(diff(log(freqsSpectrum)) / mean(diff(log(freqsSpectrum))) - 1)) < 1e-3
        spectrumFrequencyScale = 'log';
    else
        spectrumFrequencyScale = 'lin';
    end
end

%% mean and shock detection

% mean
meanWvlt = WvltComp(t, X, freqsMean, QMean, 'ZeroPadding', ZeroPadding, 'MotherWavelet', MotherWaveletMean,...
    'MeanOverFreqFunc', meanFunc);

% edge effects
[~, DeltaT] = FTpsi_DeltaT(QMean, MotherWaveletMean);
if length(ctEdgeEffectsMean) == 1
    ctEdgeEffectsMean = [ctEdgeEffectsMean, ctEdgeEffectsMean];
end
TedgeEffects = ctEdgeEffectsMean * DeltaT(freqsMean(1));
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
    figName = ['mean function : ', func2str(meanFunc), ' ; ',...
        num2str(freqsMean(1)), ' < f < ', num2str(freqsMean(end))];
    plotMeanAndShocks(t, meanWvlt, meanWvlttEdgeEffects, shockIndexes, threshold, meanTmeanWvlt, meanScale, figName);
end


%% spectrum
spectrums = WvltComp(t, X, freqsSpectrum, QSpectrum, 'ZeroPadding', ZeroPadding,...
    'MotherWavelet', MotherWaveletSpectrum, 'XindexOut', shockIndexes);
spectrums = abs(spectrums);
spectrums = transpose(spectrums);

if plotSpectrum
    plotSpectrums(freqsSpectrum, spectrums, spectrumFrequencyScale, spectrumScale);
end


end

