function [freqs, damps] = getFreqsFour(dataFileName)
%GETFREQSFOUR Summary of this function goes here
%   Detailed explanation goes here

dataFile = 'ponts marne 2\autocorr automatise\fourier\frequences.xlsx';

% lecture
T = readtable(dataFile);
freqs = T.(dataFileName).';
damps = T.([dataFileName, '_1']).';

% enlève les cases vides
damps = damps(~isnan(freqs));
freqs = freqs(~isnan(freqs));

% tri freqs croissantes
[freqs, I] = sort(freqs);
damps = damps(I);

% pourcentages
damps = damps/100;

end

