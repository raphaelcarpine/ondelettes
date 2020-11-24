function spectrums = shockDetection(t, X, freqsMean, freqsSpectrum, QMean, QSpectrum, MotherWaveletMean,...
    MotherWaveletSpectrum, ctEdgeEffectsMean,...
    meanFunc, thresholdMode, thresholdValue, maxDetectionMethod, multiSignalMean, multiSignalSpectrum, varargin)
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

%% array size

% array size
if size(t, 1) == 1 && size(t, 2) == size(X, 2)
    % ok
elseif size(t, 1) == size(X, 1) && size(t, 2) == size(X, 2)
    dt = mean(diff(t(1, :)));
    for k_t = 2:size(t, 1)
        if max(abs(t(k_t, :) - t(1, :))) > dt * 1e-3
            error(' ');
        end
    end
    t = t(1, :);
else
    error(' ');
end

% time step
if any(abs(diff(t)/mean(diff(t)) - 1) > 1e-3)
    error('non-constant time step');
end


%% mean and shock detection

% mean
meanWvltTot = nan(size(X));
for k_x = 1:size(X, 1)
    meanWvltTot(k_x, :) = WvltComp(t, X(k_x, :), freqsMean, QMean, 'MotherWavelet', MotherWaveletMean,...
        'MeanOverFreqFunc', meanFunc);
end

if multiSignalMean
    meanWvltTot = mean(meanWvltTot, 1);
end

shockIndexesTot = cell(1, size(meanWvltTot, 1));

for k_x = 1:size(meanWvltTot, 1)
    meanWvlt = meanWvltTot(k_x, :);
    
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
        if multiSignalMean
            figName = [figName, ' ; all channels'];
        else
            figName = [figName, ' ; channel ', num2str(k_x)];
        end
        plotMeanAndShocks(t, meanWvlt, meanWvlttEdgeEffects, shockIndexes,...
            threshold, meanTmeanWvlt, meanScale, figName);
        drawnow;
    end
    
    % save
    shockIndexesTot{k_x} = shockIndexes;
    
end


%% spectrum

spectrumsTot = cell(1, size(X, 1));

for k_x = 1:size(X, 1)
    if multiSignalMean
        shockIndexes = shockIndexesTot{1};
    else
        shockIndexes = shockIndexesTot{k_x};
    end
    spectrums = WvltComp(t, X(k_x, :), freqsSpectrum, QSpectrum,...
        'MotherWavelet', MotherWaveletSpectrum, 'XindexOut', shockIndexes);
    spectrums = transpose(spectrums);
    spectrumsTot{k_x} = spectrums;
end

if multiSignalMean && multiSignalSpectrum
    squaredSpectrums = zeros(size(spectrumsTot{1}));
    for k_x = 1:size(X, 1)
        squaredSpectrums = squaredSpectrums + abs(spectrumsTot{k_x}).^2;
    end
    if plotSpectrum
        figName = ['mean(abs(CWT(t=t_k, f))^2) ; all channels'];
        plotSpectrums(freqsSpectrum, squaredSpectrums, spectrumFrequencyScale, spectrumScale, figName);
    end
elseif ~multiSignalSpectrum
    for k_x = 1:size(X, 1)
        spectrums = spectrumsTot{k_x};
        if plotSpectrum
            figName = ['CWT(t=t_k, f) ; channel ', num2str(k_x)];
            plotSpectrums(freqsSpectrum, abs(spectrums), spectrumFrequencyScale, spectrumScale, figName);
        end
    end
else
    warning(' ');
end


end

