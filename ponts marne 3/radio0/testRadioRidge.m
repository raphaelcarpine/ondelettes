[T, X, dt, chDist] = getDataRadio('Mesure 5', 'mesures preliminaires');

% ~mid-span chanel
selectDist = 31.6;
Kdist = ismember(chDist, selectDist);
chDist = chDist(Kdist);
X = X(Kdist, :);

% local mean removal
Taverage = 60;
X = X - getSmoothSignal(T, X, 'gaussian', Taverage);

% ridge
freqs = linspace(1.8, 2.5, 100);
Q = 6;
MotherWavelet = 'cauchy';
signalDerivation = 2;
ct = 3;
ridgeContinuity = 'none';
slopeTimeConst = 3;
CWT = WvltComp(T, X, freqs, Q, 'MotherWavelet', MotherWavelet,...
    'DerivationOrder', signalDerivation, 'DisplayWaitBar', false);
ridge = SingleRidgeExtract(T, freqs, CWT, MotherWavelet, Q, ct, ridgeContinuity, slopeTimeConst);
Fridge = ridge.freq;
Tridge = ridge.time;
Xridge = X(ismember(T, Tridge));

%%

% plot time ridge
figure;
yyaxis left
plot(Tridge, Xridge);
ylabel('Déplacement vertical [mm]');
xlabel('Temps [s]');
yyaxis right
plot(Tridge, Fridge);
ylabel('Fréquence [Hz]');

% plot deflection ridge
figure;
plot(Xridge, Fridge);
xlabel('Déplacement vertical [mm]');
ylabel('Fréquence [Hz]');
RegressionMenu();









