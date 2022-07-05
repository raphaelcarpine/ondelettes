% folder & filename
folder = 'C:\Users\carpine\Documents\projets\ponts marne\reprise operations 2022\donnees\radio';
fileName = 'esbly0107_aprem_part3';

% read .csv
X = readtable(fullfile(folder, [fileName, '.csv']));


% data processing
chDist = {};
for k = 2:size(X, 2)
    chDist{end+1} = X.Properties.VariableNames{k};
end

X = table2array(X);
T = X(:, 1).';
X = X(:, 2:end).';

dt = (T(end) - T(1))/(length(T)-1);
if any(abs(diff(T)/dt - 1) > 1e-3)
    warning('time step pb');
end

% saving to .mat file
% save(fullfile(folder, [fileName, '.mat']), 'T', 'X', 'dt', 'chDist');


