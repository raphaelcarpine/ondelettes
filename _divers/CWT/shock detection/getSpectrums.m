function [spectrumsTot, spectrumsTot0, averageShockSpectrum0, averageSpectrum0,...
    averageUnderThresholdSpectrum0, averageAboveThresholdSpectrum0]...
    = getSpectrums(t, X, signalChannels, shockIndexes, freqs, Q, MotherWavelet, spectrumFrequencyScale, spectrumScale,...
    multiSignalMean, multiSignalSpectrum, plotShockSpectrum, varargin)
%GETSPECTRUMS Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;

plotAverageShockSpectrumDef = false;
plotAverageSpectrumDef = false;
plotUnderThresholdAverageDef = false;
plotUnderThresholdAverageIndexesDef = cell(1, size(X, 1));
plotAboveThresholdAverageDef = false;
plotAboveThresholdAverageIndexesDef = cell(1, size(X, 1));
spectrumsTot0Def = {};
averageShockSpectrum0Def = {};
averageSpectrum0Def = {};
averageUnderThresholdSpectrum0Def = {};
averageAboveThresholdSpectrum0Def = {};

addParameter(p, 'plotAverageShockSpectrum', plotAverageShockSpectrumDef);
addParameter(p, 'plotAverageSpectrum', plotAverageSpectrumDef);
addParameter(p, 'plotUnderThresholdAverage', plotUnderThresholdAverageDef);
addParameter(p, 'plotUnderThresholdAverageIndexes', plotUnderThresholdAverageIndexesDef);
addParameter(p, 'plotAboveThresholdAverage', plotAboveThresholdAverageDef);
addParameter(p, 'plotAboveThresholdAverageIndexes', plotAboveThresholdAverageIndexesDef);
addParameter(p, 'spectrumsTot0', spectrumsTot0Def);
addParameter(p, 'averageShockSpectrum0', averageShockSpectrum0Def);
addParameter(p, 'averageSpectrum0', averageSpectrum0Def);
addParameter(p, 'averageUnderThresholdSpectrum0', averageUnderThresholdSpectrum0Def);
addParameter(p, 'averageAboveThresholdSpectrum0', averageAboveThresholdSpectrum0Def);

parse(p, varargin{:});

plotAverageShockSpectrum = p.Results.plotAverageShockSpectrum;
plotAverageSpectrum = p.Results.plotAverageSpectrum;
plotUnderThresholdAverage = p.Results.plotUnderThresholdAverage;
plotUnderThresholdAverageIndexes = p.Results.plotUnderThresholdAverageIndexes;
plotAboveThresholdAverage = p.Results.plotAboveThresholdAverage;
plotAboveThresholdAverageIndexes = p.Results.plotAboveThresholdAverageIndexes;
spectrumsTot0 = p.Results.spectrumsTot0;
averageShockSpectrum0 = p.Results.averageShockSpectrum0;
averageSpectrum0 = p.Results.averageSpectrum0;
averageUnderThresholdSpectrum0 = p.Results.averageUnderThresholdSpectrum0;
averageAboveThresholdSpectrum0 = p.Results.averageAboveThresholdSpectrum0;

%%

if multiSignalMean
    shockIndexes1 = shockIndexes;
    shockIndexes = cell(1, size(X, 1));
    shockIndexes(:) = shockIndexes1;
    
    plotUnderThresholdAverageIndexes1 = plotUnderThresholdAverageIndexes;
    plotUnderThresholdAverageIndexes = cell(1, size(X, 1));
    plotUnderThresholdAverageIndexes(:) = plotUnderThresholdAverageIndexes1;
    
    plotAboveThresholdAverageIndexes1 = plotAboveThresholdAverageIndexes;
    plotAboveThresholdAverageIndexes = cell(1, size(X, 1));
    plotAboveThresholdAverageIndexes(:) = plotAboveThresholdAverageIndexes1;
end


%% spectrums computing

% computing CWT
if isempty(spectrumsTot0) || isempty(averageShockSpectrum0) || isempty(averageSpectrum0) ||...
        isempty(averageUnderThresholdSpectrum0) || isempty(averageAboveThresholdSpectrum0)
    spectrumsTot0 = cell(1, size(X, 1)); % shock spectrums
    averageShockSpectrum0 = cell(1, size(X, 1)); % average shock spectrum
    averageSpectrum0 = cell(1, size(X, 1)); % average spectrum
    averageUnderThresholdSpectrum0 = cell(1, size(X, 1)); % average under threshold spectrum
    averageAboveThresholdSpectrum0 = cell(1, size(X, 1)); % average above threshold spectrum
    
    for k_ch = 1:size(X, 1)
        MeanSquareOverXelements = {shockIndexes{k_ch},... % average shock spectrum
            true(1, size(X, 2)),... % average spectrum
            plotUnderThresholdAverageIndexes{k_ch},... % average under threshold spectrum
            plotAboveThresholdAverageIndexes{k_ch}}; % average above threshold spectrum
        spectrums = WvltComp(t, X(k_ch, :), freqs, Q, 'MotherWavelet', MotherWavelet, 'XindexOut',...
            shockIndexes{k_ch}, 'MeanSquareOverXelements', MeanSquareOverXelements);
        spectrums = transpose(spectrums);
        
        averageShockSpectrum0{k_ch} = spectrums(end-3, :); % average shock spectrum
        averageSpectrum0{k_ch} = spectrums(end-2, :); % average spectrum
        averageUnderThresholdSpectrum0{k_ch} = spectrums(end-1, :); % average under threshold spectrum
        averageAboveThresholdSpectrum0{k_ch} = spectrums(end, :); % average above threshold spectrum
        
        spectrums = spectrums(1:end-4, :);
        spectrums = abs(spectrums).^2;
        spectrumsTot0{k_ch} = spectrums;
    end
end


spectrumsTot = spectrumsTot0;
averageShockSpectrum = averageShockSpectrum0;
averageSpectrum = averageSpectrum0;
averageUnderThresholdSpectrum = averageUnderThresholdSpectrum0;
averageAboveThresholdSpectrum = averageAboveThresholdSpectrum0;

if multiSignalSpectrum % averaging over multiple channels
    spectrumsTot = {sum(cat(3, spectrumsTot{:}), 3) / length(spectrumsTot)};
    averageShockSpectrum = {sum(cat(3, averageShockSpectrum{:}), 3) / length(averageShockSpectrum)};
    averageSpectrum = {sum(cat(3, averageSpectrum{:}), 3) / length(averageSpectrum)};
    averageUnderThresholdSpectrum = {sum(cat(3, averageUnderThresholdSpectrum{:}), 3) / length(averageUnderThresholdSpectrum)};
    averageAboveThresholdSpectrum = {sum(cat(3, averageAboveThresholdSpectrum{:}), 3) / length(averageAboveThresholdSpectrum)};
end



%% plot

% plotShockSpectrum
% plotAverageShockSpectrum
% plotAverageSpectrum
% plotUnderThresholdAverage
% plotAboveThresholdAverage

deltaFwvlt = @(f) f/(2*Q);

% all shocks plot
if plotShockSpectrum
    if multiSignalSpectrum
        figName = '|CWT(t=t_k, f)|² ; all selected channels';
        plotSpectrums(freqs, spectrumsTot{1}, spectrumFrequencyScale, spectrumScale, deltaFwvlt, figName);
    else
        for k_ch = 1:length(spectrumsTot)
            figName = sprintf('|CWT(t=t_k, f)|² ; channel %d', signalChannels(k_ch));
            plotSpectrums(freqs, spectrumsTot{k_ch}, spectrumFrequencyScale, spectrumScale, deltaFwvlt, figName);
        end
    end
end

% averages plot
if plotAverageShockSpectrum || plotAverageSpectrum || plotUnderThresholdAverage || plotAboveThresholdAverage
    % spectrums array
    spectrumsArray = cell(size(averageShockSpectrum));
    for k_ch = 1:length(averageShockSpectrum)
        spectrumsArray{k_ch} = [averageShockSpectrum{k_ch};
            averageSpectrum{k_ch};
            averageUnderThresholdSpectrum{k_ch};
            averageAboveThresholdSpectrum{k_ch}];
    end
    
    % plot
    if multiSignalSpectrum
        figName = 'mean(|CWT(~, f)|²) ; all selected channels';
        plotSpectrums(freqs, spectrumsArray{1}, spectrumFrequencyScale, spectrumScale, deltaFwvlt, figName,...
            plotAverageShockSpectrum, plotAverageSpectrum, plotUnderThresholdAverage, plotAboveThresholdAverage);
    else
        for k_ch = 1:length(spectrumsTot)
            figName = sprintf('mean(|CWT(~, f)|²) ; channel %d', signalChannels(k_ch));
            plotSpectrums(freqs, spectrumsArray{k_ch}, spectrumFrequencyScale, spectrumScale, deltaFwvlt, figName,...
                plotAverageShockSpectrum, plotAverageSpectrum, plotUnderThresholdAverage, plotAboveThresholdAverage);
        end
    end
end


end

