
filePath = getResultsFile(12);

load(filePath);

[~, fileName] = fileparts(filePath);
disp(fileName);

var_temp = @(t) (deltaT/2) * sin(2*pi*t/periodeTemp);
Temp = var_temp(T);

%% options

methodeFreq = 'CWT'; % 'CWT', 'hilbert'
intAccPos = 0; % integration acceleration pour position

% hilbert
fc1 = 0.5;
fc2 = 5;

% CWT
fmin = 1;
fmax = 3;
Q = 2;
MotherWavelet = 'morlet';
ct = 3;
ridgeContinuity = 0;

plotFA = 1;
plotK = 0;
plotDiscardedFA = 1;
plotRegLin = 1;

referenceQtyArray = {'ampl', 'pos', 'temp'};
averageScale = 'lin';
averageIncrementArray = [0.001, 0.0005, 0.5];

filtrageSignal = 0;
filtragePos = 0;
discardThresholdAmpl = 0;

averageTimeIntervals = 25;
Kmin = 10;

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

% filtrage signal acc
if filtrageSignal
    figure;
    plot(T, Acapt);
    Acapt = butterworthFilter(T, Acapt, 4, 'low', 5);
    hold on
    plot(T, Acapt);
end

% filtrage pos
if filtragePos
    figure;
    plot(T, Ycapt);
    Ycapt = butterworthFilter(T, Ycapt, 1, 'low', 5);
%     Ycapt = getSmoothSignal(T, Ycapt, 'gaussian', 1.5);
    hold on
    plot(T, Ycapt);
end

%% position, en intégrant 2 fois

if intAccPos
    windowType = 'gaussian'; % fenpetre de lissage
    Tmean = 50; % largeur fenêtre de lissage
    Amean = Acapt - getSmoothSignal(T, Acapt, windowType, Tmean);
    IA = mean(diff(T)) * cumsum(Amean); % integration
    IAmean = IA - getSmoothSignal(T, IA, windowType, Tmean);
    IIA = mean(diff(T)) * cumsum(IAmean); % integration
    Ycapt = IIA;
end

%% analyse CWT ou hilbert

% recherche si ridges déjà calculés
if strcmp(methodeFreq, 'CWT')
    ridgeFolderPath = 'C:\Users\carpine\Documents\MATLAB\ondelettes\simulations el finis nonlin\data ridges';
    continuouStr = '_continuous';
    ridgeFolder = sprintf('ridges_fmin%g_fmax%g_Q%g_%s%s', [fmin, fmax, Q,...
        convertCharsToStrings(MotherWavelet), convertCharsToStrings(continuouStr(1:end*ridgeContinuity))]);
    ridgeFileCompletePath = fullfile(ridgeFolderPath, ridgeFolder, fileName);
    if ~isfile(ridgeFileCompletePath)
        error('file not found');
    end
end


if strcmp(methodeFreq, 'CWT')
    freqs = linspace(fmin, fmax, 100);
    
    CWT = WvltComp(T, Acapt, freqs, Q, 'MotherWavelet', MotherWavelet);
    [~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
    
    indexes_ridge = T >= T(1)+ct*DeltaT(fmin) & T <= T(end)-ct*DeltaT(fmin);
    Tridge = T(indexes_ridge);
    CWTridge = [nan(1, length(Tridge)); CWT(:, indexes_ridge); nan(1, length(Tridge))];
    CWTridge0 = [zeros(1, length(Tridge)); CWT(:, indexes_ridge); zeros(1, length(Tridge))];
    freqs = [nan, freqs, nan];
    if ridgeContinuity
        Fridge = nan(1, length(Tridge));
        [~, Fridge(1)] = max(abs(CWTridge0(:, 1)));
        localMax = abs(CWTridge0(1:end-1, :)) < abs(CWTridge0(2:end, :));
        localMax = localMax(1:end-1, :) & ~localMax(2:end, :);
        for kt = 2:length(Tridge)
            localMaxFreq = find(localMax(:, kt)) + 1;
            [~, closestLocalMax] = min(abs(localMaxFreq - Fridge(kt-1)));
            Fridge(kt) = localMaxFreq(closestLocalMax);
        end
    else
        [~, Fridge] = max(abs(CWTridge), [], 1);
    end
    [Fridge, Aridge] = localMax3Points(freqs([Fridge-1; Fridge; Fridge+1]),...
        CWTridge([Fridge-1; Fridge; Fridge+1] + [1;1;1] * (0:size(CWTridge, 2)-1)*size(CWTridge, 1)));
    Aridge = abs(Aridge);
    Yridge = Ycapt(indexes_ridge);
    Tempridge = Temp(indexes_ridge);
elseif strcmp(methodeFreq, 'hilbert')
    if fc1 > 0
        Acapt = butterworthFilter(T, Acapt, fc1, 'high', 5);
    end
    if fc2 < inf
        Acapt = butterworthFilter(T, Acapt, fc2, 'low', 5);
    end
    Aridge = hilbert(Acapt);
    dt = (T(end)-T(1))/(length(T)-1);
    Fridge = angle(Aridge(2:end) ./ Aridge(1:end-1))/(2*pi*dt);
    Fridge = [Fridge(1), (Fridge(1:end-1) + Fridge(2:end))/2, Fridge(end)];
    Aridge = abs(Aridge);
    Tridge = T;
    Yridge = Ycapt;
    Tempridge = Temp;
else
    error(' ');
end

%% post traitement

% supression des valeurs nan
notNaN = ~isnan(Tridge) & ~isnan(Fridge) & ~isnan(Aridge) & ~isnan(Yridge) & ~isnan(Tempridge);
Tridge = Tridge(notNaN);
Fridge = Fridge(notNaN);
Aridge = Aridge(notNaN);
Yridge = Yridge(notNaN);
Tempridge = Tempridge(notNaN);

% suppression ridge amplitude trop faible
thresholdAmpl = 0.01;
if discardThresholdAmpl
    amplOK = Aridge >= thresholdAmpl;
    Fridge = Fridge(amplOK);
    Aridge = Aridge(amplOK);
    Yridge = Yridge(amplOK);
    Tempridge = Tempridge(amplOK);
end

%% reg lin multi dimension

coeffsRegLin = [ones(size(Fridge)); Aridge - mean(Aridge); Yridge - mean(Yridge); Tempridge - mean(Tempridge)].' \ Fridge.';

fprintf('Reg. lin. multi-param :\nf = %.4f + %.3g*(A-A0) + %.3g*(Y-Y0) + %.3g*(T-T0)\n', coeffsRegLin);
fprintf('%.3g*std(A) = %.3g Hz\n', abs(coeffsRegLin(2)) * [1, std(Aridge)]);
fprintf('%.3g*std(Y) = %.3g Hz\n', abs(coeffsRegLin(3)) * [1, std(Yridge)]);
fprintf('%.3g*std(T) = %.3g Hz\n', abs(coeffsRegLin(4)) * [1, std(Tempridge)]);


%% averaging

for kqty = 1:length(referenceQtyArray)
    % choix quantité de reference
    if strcmp(referenceQtyArray{kqty}, 'ampl')
        Qtyridge = Aridge;
        referenceQty = 'Amplitude [m/s²]';
    elseif strcmp(referenceQtyArray{kqty}, 'pos')
        Qtyridge = Yridge;
        referenceQty = 'Position [m]';
    elseif strcmp(referenceQtyArray{kqty}, 'temp')
        Qtyridge = Tempridge;
        referenceQty = 'Temperature [°C]';
    end
    averageIncrement = averageIncrementArray(kqty);
    
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
            plot(Xax, coeffsRegLin(1) + coeffsRegLin(1+kqty)*(Xax - mean(Qtyridge)), 'r--');
        end
    end

end

%% reg lin pos par heure

coeffsPos = nan(1, 24);

for kh = 1:24
    indexesHeure = (Tridge >= (kh-1)*3600) & (Tridge < kh*3600);
    coeffsRegLin = [ones(size(Fridge(indexesHeure))); Yridge(indexesHeure) - mean(Yridge(indexesHeure))].' \ Fridge(indexesHeure).';
    coeffsPos(kh) = coeffsRegLin(2);
end


fprintf('\nReg. lin. pos :\ncoeff = %.1f +- %.1f Hz/m\n', [mean(coeffsPos), std(coeffsPos)/sqrt(24)]);
