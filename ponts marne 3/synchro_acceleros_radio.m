clear all

[filePath, fileName, filePathRadio, fileNameRadio] = choixData(0);

% acceleros
load(filePath);
A = X.';
T = T.';
[A, T] = removeRedundantData(A, T);
[A, T, DeltaT] = removeNanSignal(A, T);

k1 = find(contains(channelNames, '29280:ch3'));
k2 = find(contains(channelNames, '40199:ch3'));
A = A(k1, :) + A(k2, :);


% radio
load(filePathRadio);

Xradio(isnan(Xradio)) = 0;
X = Xradio(2, :);
% X = diff(diff(Xradio(2, :)));


%% xcorr

[r, lags] = xcorr(A, X, length(T)-length(Tradio));

figure('Name', fileName);
lags0 = find(lags >= 0);
plt = plot(lags(lags0)/128, r(lags>=0));
xlim([0, lags(lags0(end))]/128);
xlabel('Délai [s]');
ylabel('Intercorrélation (accéléromètres, interféromètre)');
plt.DataTipTemplate.DataTipRows = plt.DataTipTemplate.DataTipRows(1);
plt.DataTipTemplate.DataTipRows.Label = 'k =';
plt.DataTipTemplate.DataTipRows.Value = lags(lags0);
plt.DataTipTemplate.DataTipRows.Format = 'auto';

% xline((seconds(startDateRadio - startDate) - 0.5 - DeltaT), 'LineWidth', 1);
% xline((seconds(startDateRadio - startDate) + 1 - DeltaT), 'LineWidth', 1);
xline((seconds(startDateRadio - startDate) + 0.5 - DeltaT), 'k', 'LineWidth', 1,...
    'Label', 'HEURE PC', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');

% changis : 467
% esbly matin : 476
% esbly aprem : 1810













