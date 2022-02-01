clear all

dataFolder = 'ponts marne 2\erreur amort\save';

Nsep = 20;


[~, dataFileName] = fileparts(choixData(6));

load(fullfile(dataFolder, ['modes_', dataFileName, '_Nsep', num2str(Nsep)]));

Freqs0 =      [2.08 2.18 2.84 2.94 5.47 6.52 7.85 15.04 16.74 21.13];
Damps0 = 0.01*[0.9  0.9  0.8  1.0  1.8  3.7  5.2  2.4   1.0   1.0  ];
Poles0 = 2*pi * Freqs0 .*( - Damps0 + 1i * sqrt(1 - Damps0.^2));
Topt = -1.4326./real(Poles0);

Meff0 = 1;

%%

Poles = 2*pi * Freqs .*( - Damps + 1i * sqrt(1 - Damps.^2));

meanFreqs = mean(Freqs, 2, 'omitnan');
stdFreqs = std(Freqs, 0, 2, 'omitnan');
meanDamps = 100*mean(Damps, 2, 'omitnan');
stdDamps = 100*std(Damps, 0, 2, 'omitnan');

Tinit = Tinit(:, 1);

alpha = sqrt((-4./real(Poles0)+2./imag(Poles0))/Tsep).';
Meff = (0.01*stdDamps.' ./ Damps0 .* exp(real(Poles0).*Tinit.') ./ alpha.').';
alphaI = alpha .* exp(-real(Poles0).'.*Tinit);
alphaF = alpha .* exp(-real(Poles0).'.*(Tinit+Topt.'));

expectedStdDamps = 100*Damps0.' * Meff0 .* exp(-real(Poles0).'.*Tinit) .* alpha;
expectedMeanDamps = 100*Damps0.';
expectedCVDamps = Meff0 .* exp(-real(Poles0).'.*Tinit) .* alpha;

T = table(meanFreqs, stdFreqs, meanDamps, expectedMeanDamps, stdDamps, expectedStdDamps, alpha, alphaI, alphaF, Meff);
% T

modes = (1:10).';

T2 = table(modes, alphaF, meanFreqs, meanDamps, expectedStdDamps);
% T2

T3 = table(modes, alphaF, meanFreqs, meanDamps, stdDamps, Meff);
T3


