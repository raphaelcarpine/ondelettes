function getSpectrums(t, X, shockIndexes, freqs, Q, MotherWavelet, spectrumFrequencyScale, spectrumScale,...
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

addParameter(p, 'plotAverageShockSpectrum', plotAverageShockSpectrumDef);
addParameter(p, 'plotAverageSpectrum', plotAverageSpectrumDef);
addParameter(p, 'plotUnderThresholdAverage', plotUnderThresholdAverageDef);
addParameter(p, 'plotUnderThresholdAverageIndexes', plotUnderThresholdAverageIndexesDef);
addParameter(p, 'plotAboveThresholdAverage', plotAboveThresholdAverageDef);
addParameter(p, 'plotAboveThresholdAverageIndexes', plotAboveThresholdAverageIndexesDef);

parse(p, varargin{:});

plotAverageShockSpectrum = p.Results.plotAverageShockSpectrum;
plotAverageSpectrum = p.Results.plotAverageSpectrum;
plotUnderThresholdAverage = p.Results.plotUnderThresholdAverage;
plotUnderThresholdAverageIndexes = p.Results.plotUnderThresholdAverageIndexes;
plotAboveThresholdAverage = p.Results.plotAboveThresholdAverage;
plotAboveThresholdAverageIndexes = p.Results.plotAboveThresholdAverageIndexes;

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


spectrumsTot = cell(1, size(X, 1)); % shock spectrums
averageShockSpectrum = cell(1, size(X, 1)); % average shock spectrum
averageSpectrum = cell(1, size(X, 1)); % average spectrum
averageUnderThresholdSpectrum = cell(1, size(X, 1)); % average under threshold spectrum
averageAboveThresholdSpectrum = cell(1, size(X, 1)); % average above threshold spectrum

for k_ch = 1:size(X, 1)
    MeanSquareOverXelements = {shockIndexes{k_ch},... % average shock spectrum 
        true(1, size(X, 2)),... % average spectrum
        plotUnderThresholdAverageIndexes{k_ch},... % average under threshold spectrum
        plotAboveThresholdAverageIndexes{k_ch}}; % average above threshold spectrum
    spectrums = WvltComp(t, X(k_ch, :), freqs, Q, 'MotherWavelet', MotherWavelet, 'XindexOut',...
        shockIndexes{k_ch}, 'MeanSquareOverXelements', MeanSquareOverXelements);
    spectrums = transpose(spectrums);
    
    averageShockSpectrum{k_ch} = spectrums(end-3, :); % average shock spectrum
    averageSpectrum{k_ch} = spectrums(end-2, :); % average spectrum
    averageUnderThresholdSpectrum{k_ch} = spectrums(end-1, :); % average under threshold spectrum
    averageAboveThresholdSpectrum{k_ch} = spectrums(end, :); % average above threshold spectrum
    
    spectrums = spectrums(1:end-4, :);
    spectrums = abs(spectrums).^2;
    spectrumsTot{k_ch} = spectrums;
end

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

% all shocks plot
if plotShockSpectrum
    if multiSignalSpectrum
        figName = '|CWT(t=t_k, f)|² ; all channels';
        plotSpectrums(freqs, spectrumsTot{1}, spectrumFrequencyScale, spectrumScale, figName);
    else
        for k_ch = 1:length(spectrumsTot)
            figName = sprintf('|CWT(t=t_k, f)|² ; channel %d', k_ch);
            plotSpectrums(freqs, spectrumsTot{k_ch}, spectrumFrequencyScale, spectrumScale, figName);
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
        figName = 'mean(|CWT(~, f)|²) ; all channels';
        plotSpectrums(freqs, spectrumsArray{1}, spectrumFrequencyScale, spectrumScale, figName,...
            plotAverageShockSpectrum, plotAverageSpectrum, plotUnderThresholdAverage, plotAboveThresholdAverage);
    else
        for k_ch = 1:length(spectrumsTot)
            figName = sprintf('mean(|CWT(~, f)|²) ; channel %d', k_ch);
            plotSpectrums(freqs, spectrumsArray{k_ch}, spectrumFrequencyScale, spectrumScale, figName,...
                plotAverageShockSpectrum, plotAverageSpectrum, plotUnderThresholdAverage, plotAboveThresholdAverage);
        end
    end
end


end

