
filePath = getResultsFile(-1);

load(filePath);

Temp = var_temp(T);

%% options

plotFA = 1;
plotK = 0;
plotDiscardedFA = 1;
plotRegLin = 1;

referenceQty = 'pos'; % 'ampl', 'pos', 'temp'
averageScale = 'lin';
averageTimeIntervals = 25;
averageIncrement = .0001;

Kmin = 3;

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
Q = 6;
MotherWavelet = 'morlet';
ct = 3;

freqs = linspace(fmin, fmax, 300);

CWT = WvltComp(T, Acapt, freqs, Q, 'MotherWavelet', MotherWavelet);
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
Tempridge = Temp(indexes_ridge);

% choix quantité de reference
if strcmp(referenceQty, 'ampl')
    Qtyridge = Aridge;
    referenceQty = 'Amplitude [m/s²]';
elseif strcmp(referenceQty, 'pos')
    Qtyridge = Yridge;
    referenceQty = 'Position [m]';
elseif strcmp(referenceQty, 'temp')
    Qtyridge = Tempridge;
    referenceQty = 'Temperature [°C]';
end

%% reg lin

Qtyridge2 = Qtyridge(~isnan(Qtyridge) & ~isnan(Fridge));
Fridge2 = Fridge(~isnan(Qtyridge) & ~isnan(Fridge));

coeffsRegLin = [ones(size(Qtyridge2)); Qtyridge2].' \ Fridge2.';


%% averaging

[X, Y, stdY, K, Xlims] = averagingScatter3(Qtyridge, Fridge, averageIncrement, averageScale, Tridge, averageTimeIntervals);

Y0 = Y;
stdY0 = stdY;
Y(K < Kmin) = nan;
stdY(K < Kmin) = nan;

if plotK
    figure('Name', '');
    if strcmp(averageScale, 'lin')
        bar(X, K);
    else
        bar(log10(X), K);
        set(gca,'Xtick',-5:1);
        set(gca,'Xticklabel',10.^get(gca,'Xtick'));
    end
    set(gca, 'XScale', averageScale);
    xlabel(referenceQty);
end

if plotFA
    figure('Name', '');
    if plotDiscardedFA
        err0 = 1.96*stdY0./sqrt(K);
        l0 = errorbar(X, Y0, err0);
        l0.Color = [1 1 1] - 0.5*([1 1 1]-l0.Color);
        hold on
        ax = gca;
        ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
    end
    err = 1.96*stdY./sqrt(K);
    errorbar(X, Y, err);
    set(gca, 'XScale', averageScale);
    xlabel(referenceQty);
    ylabel('Frequency [Hz]');
    if plotRegLin
        hold on
        Xax = get(gca, 'XLim');
        Xax = linspace(Xax(1), Xax(2));
        plot(Xax, coeffsRegLin(1) + coeffsRegLin(2)*Xax, 'r--');
    end
end

