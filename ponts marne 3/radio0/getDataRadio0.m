function [T, X, dt, chDist] = getDataRadio(fileName, subFolder)
%GETDATARADIO Summary of this function goes here
%   Detailed explanation goes here

% folder & filename
folder = 'C:\Users\carpine\Documents\projets\ponts marne\reprise opérations 2022\mesures radio\data';
if nargin >= 2
    folder = fullfile(folder, subFolder);
end

[~, fileName] = fileparts(fileName);

% load .mat file if exists
try
    load(fullfile(folder, [fileName, '.mat']));
    return
catch
end

% read .xlsx if not
X = readmatrix(fullfile(folder, [fileName, '.xlsx']));

% data processing
chDist = X(1, 2:end);
k = 1;
while isnan(X(k, 1))
    k = k + 1;
end
T = X(k:end, 1).';
X = X(k:end, 2:end).';

dt = (T(end) - T(1))/(length(T)-1);
if any(abs(diff(T)/dt - 1) > 1e-3)
    warning('time step pb');
end

% saving to .mat file
save(fullfile(folder, [fileName, '.mat']), 'T', 'X', 'dt', 'chDist');

end

