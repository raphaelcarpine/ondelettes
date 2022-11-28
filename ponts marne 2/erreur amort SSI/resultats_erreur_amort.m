clear all

dataFolder = 'ponts marne 2\erreur amort SSI\save';

Nsep = 20;


[~, dataFileName] = fileparts(choixData(6));

load(fullfile(dataFolder, ['modes_', dataFileName, '_Nsep', num2str(Nsep)]));

Freqs0 =      [2.08 2.18 2.84 2.94 5.47 6.52 7.85 15.04 16.74 21.13];
Damps0 = 0.01*[0.9  0.9  0.8  1.0  1.8  3.7  5.2  2.4   1.0   1.0  ];
Poles0 = 2*pi * Freqs0 .*( - Damps0 + 1i * sqrt(1 - Damps0.^2));


%%

Poles = 2*pi * Freqs .*( - Damps + 1i * sqrt(1 - Damps.^2));

meanFreqs = mean(Freqs, 2, 'omitnan');
stdFreqs = std(Freqs, 0, 2, 'omitnan');
meanDamps = 100*mean(Damps, 2, 'omitnan');
stdDamps = 100*std(Damps, 0, 2, 'omitnan');

T = table(meanFreqs, stdFreqs, meanDamps, stdDamps);
Ts


