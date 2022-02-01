
filePath = getResultsFile();

load(filePath);

%% calcul, mise en forme matrices

% interpolation etc

getYtot = @(Y) [zeros(1, size(Y, 2)); Y; zeros(1, size(Y, 2))];

% DDL 1, 2, end-1, end
Ytot = getYtot(Y);
Vtot = getYtot(V);
Atot = getYtot(A);

% capteurs
pos_capteurs = L/2;
Ycapt = getYcapt2(Ytot, pos_capteurs, dx);
Vcapt = getYcapt2(Vtot, pos_capteurs, dx);
Acapt = getYcapt2(Atot, pos_capteurs, dx);

%% analyse CWT

fmin = 1;
fmax = 3;
Q = 4;
MotherWavelet = 'morlet';
ct = 3;

freqs = linspace(fmin, fmax, 300);

CWT = WvltComp(T, Ycapt, freqs, Q, 'MotherWavelet', MotherWavelet);
[~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);

indexes_ridge = T >= T(1)+ct*DeltaT(fmin) & T <= T(end)-ct*DeltaT(fmin);
Tridge = T(indexes_ridge);
CWTridge = [nan(1, length(Tridge)); CWT(:, indexes_ridge); nan(1, length(Tridge))];
freqs = [nan, freqs, nan];
[~, Fridge] = max(abs(CWTridge), [], 1);
[Fridge, Aridge] = localMax3Points(freqs([Fridge-1; Fridge; Fridge+1]),...
    CWTridge([Fridge-1; Fridge; Fridge+1] + [1;1;1] * (0:size(CWTridge, 2)-1)*size(CWTridge, 1)));
Aridge = abs(Aridge);
Yridge = Ycapt(indexes_ridge);

[X, Y, stdY, K, Xlims] = averagingScatter3(Aridge, Fridge, 0.01, averageScale, Tridge, 100);

end