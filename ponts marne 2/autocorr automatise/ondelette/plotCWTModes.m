clear all

dataFolder = 'ponts marne 2\autocorr automatise\ondelette\save';

[~, dataFileName] = fileparts(choixData(0));

load(fullfile(dataFolder, ['modes_', dataFileName]));

try
    Freqs = Freqs1;
    Damps  = Damps1;
    Shapes = Shapes1;
catch
end

dimensionsShapes2;

for kf = 1:length(Freqs)
    figName = sprintf('f = %.2fHz, z = %.2f%%; In = %.2f%%', [Freqs(kf), 100*Damps(kf), 100*nonPropIndex(Shapes(:, kf))]);
    figName = [dataFileName, '; ', figName];
    if PbCalculRidge(kf)
        figName = [figName, ' (pb extractio nridge)'];
    end
    shapePlotBridge(real(Shapes(:, kf)), figName);
%     shapePlotBridgeAnim(Shapes(:, kf), figName);
end