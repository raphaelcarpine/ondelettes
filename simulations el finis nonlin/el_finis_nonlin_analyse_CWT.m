clear all

%% options

computeAll = 0;

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
ridgeContinuity = 'none'; % 'none', 'simple', 'slope0.1', 'slope0.3', 'slope1'
signalDerivation = 0;

% linear regressions
filtrageSignal = 0;
filtragePos = 0;
discardThresholdAmpl = 0;
thresholdAmpl = 0.01;

regLiMultiVarsNames = {'A', 'Y', 'T'};
regLiMultiVars = [1 1 1]; %ampl, pos, temp
regLiMultiVars = [true, logical(regLiMultiVars)];

% averaging
averageTimeIntervals = 25;
Kmin = 10;

plotFA = 0;
plotK = 0;
plotDiscardedFA = 1;
plotRegLin = 1;

referenceQtyArray = {'ampl', 'pos', 'temp'};
averageScale = 'lin';
averageIncrementArray = [0.002*(2*pi*2.08)^signalDerivation, 0.0005, 0.5];

%% folders
if exist('D:\simulations elements finis non lin\data', 'dir')
    ridgeFolderPath = 'D:\simulations elements finis non lin\data\data ridges';
elseif exist('C:\Users\raphael\Documents\resultats simul diff finies', 'dir')
    ridgeFolderPath = 'C:\Users\raphael\Documents\resultats simul diff finies\data ridges';
else
    error(' ');
end

if ~computeAll
    signalDerivationStrs = ["_2integration", "_integration", "", "_derivation", "_2derivation"];
    ridgeFolder = sprintf('ridges_fmin%g_fmax%g_Q%g_%s%s_%s', [fmin, fmax, Q,...
        convertCharsToStrings(MotherWavelet), signalDerivationStrs(signalDerivation+3), convertCharsToStrings(ridgeContinuity)]);
    ridgeFolders = {ridgeFolder};
else
    folders = dir(ridgeFolderPath);
    ridgeFolders = {folders.name};
    ridgeFolders = ridgeFolders([folders.isdir]);
    ridgeFolders = ridgeFolders(3:end);
end

%%

for kfolder = 1:length(ridgeFolders)
    ridgeFolder = ridgeFolders{kfolder};

    folderInfo = strsplit(ridgeFolder, '_');
    Q = str2double(folderInfo{4}(2:end));
    MotherWavelet = folderInfo{5};


    Ksimu = [1:36, 101, 106].';
    F0 = nan(size(Ksimu));
    F0_err = nan(size(Ksimu));
    coeffAmpl = nan(size(Ksimu));
    coeffAmpl_err = nan(size(Ksimu));
    coeffPos = nan(size(Ksimu));
    coeffPos_err = nan(size(Ksimu));
    coeffTemp = nan(size(Ksimu));
    coeffTemp_err = nan(size(Ksimu));

    for ks = 1:length(Ksimu)
        ksim = Ksimu(ks);
        filePath = getResultsFile(ksim);

        load(filePath, 'A', 'Y', 'L', 'T', 'dx', 'deltaT', 'periodeTemp');

        [~, fileName] = fileparts(filePath);
        disp(fileName);

        var_temp = @(t) (deltaT/2) * sin(2*pi*t/periodeTemp);
        Temp = var_temp(T);


        %% calcul, mise en forme matrices

        % interpolation etc

        getYtot = @(Y) [zeros(1, size(Y, 2)); Y; zeros(1, size(Y, 2))];

        % DDL 1, 2, end-1, end
        Ytot = getYtot(Y);
        %     Vtot = getYtot(V);
        Atot = getYtot(A);

        % capteurs
        pos_capteurs = L/2;
        Ycapt = getYcapt2(Ytot, pos_capteurs, dx);
        %     Vcapt = getYcapt2(Vtot, pos_capteurs, dx);
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


        if strcmp(methodeFreq, 'CWT')
            % load ridge file
            load(fullfile(ridgeFolderPath, ridgeFolder, [fileName, '.mat']));

            % edge effects
            [~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
            indexes_ridge = T >= T(1)+ct*DeltaT(fmin) & T <= T(end)-ct*DeltaT(fmin);
            Tridge = T(indexes_ridge);
            Fridge = Fridge(indexes_ridge);
            Aridge = abs(Aridge(indexes_ridge));
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
        if discardThresholdAmpl
            amplOK = Aridge >= thresholdAmpl;
            Tridge = Tridge(amplOK);
            Fridge = Fridge(amplOK);
            Aridge = Aridge(amplOK);
            Yridge = Yridge(amplOK);
            Tempridge = Tempridge(amplOK);
        end

        %% reg lin multi dimension

        regLiMultiVarsVals = [ones(size(Fridge)); Aridge; Yridge; Tempridge];
        regLiMultiVarsVals = regLiMultiVarsVals(regLiMultiVars, :);
        RLMVmult = ones(1, size(regLiMultiVarsVals, 1)); % mean selection, to avoid removing mean of 1
        RLMVmult(1) = 0;
        RLMVmult(1:end) = 0; % no mean removal

        [coeffsRegLin0, coeffsRegLinConfIntervals0] =...
            regress(Fridge.', regLiMultiVarsVals.' - RLMVmult.*mean(regLiMultiVarsVals.'));
        coeffsRegLin = zeros(length(regLiMultiVars), 1);
        coeffsRegLinConfIntervals = zeros(length(regLiMultiVars), 2);
        coeffsRegLin(regLiMultiVars) = coeffsRegLin0;
        coeffsRegLinConfIntervals(regLiMultiVars, :) = coeffsRegLinConfIntervals0;
        coeffsRegLinConfIntervals = (coeffsRegLinConfIntervals(:, 2) - coeffsRegLinConfIntervals(:, 1))/2;

        % fprintf('Reg. lin. multi-param :\nf = %.4f + %.3g*(A-A0) + %.3g*(Y-Y0) + %.3g*(T-T0)\n', coeffsRegLin);
        fprintf('Reg. lin. multi-param :\nf = (%.4f+-%.4f) + (%.3g+-%.3g)*(A-A0) + (%.3g+-%.3g)*(Y-Y0) + (%.3g+-%.3g)*(T-T0)\n',...
            [coeffsRegLin.'; coeffsRegLinConfIntervals.']);
        fprintf('%.3g*std(A) = %.3g Hz\n', abs(coeffsRegLin(2)) * [1, std(Aridge)]);
        fprintf('%.3g*std(Y) = %.3g Hz\n', abs(coeffsRegLin(3)) * [1, std(Yridge)]);
        fprintf('%.3g*std(T) = %.3g Hz\n\n', abs(coeffsRegLin(4)) * [1, std(Tempridge)]);

        % save
        F0(ks) = coeffsRegLin(1);
        coeffAmpl(ks) = coeffsRegLin(2);
        coeffPos(ks) = coeffsRegLin(3);
        coeffTemp(ks) = coeffsRegLin(4);

        %% reg lin multi dimension par minute de chaque heure
        
        Ndecoup = 20; % decoupe heure

        coeffsRegLinMin = zeros(4, Ndecoup);

        for kmin = 1:Ndecoup
            indexesMin = mod(floor(Tridge*Ndecoup/3600), Ndecoup) == kmin-1;
            coeffsRegLinMin(regLiMultiVars, kmin) = ...
                (regLiMultiVarsVals(:, indexesMin).' - RLMVmult.*mean(regLiMultiVarsVals(:, indexesMin).'))\ Fridge(indexesMin).';
        end

        fprintf('Reg. lin. multi-param / min :\nf = (%.4f+-%.4f) + (%.3g+-%.3g)*(A-A0) + (%.3g+-%.3g)*(Y-Y0) + (%.3g+-%.3g)*(T-T0)\n',...
            [mean(coeffsRegLinMin.'); 1.96*std(coeffsRegLinMin.')/sqrt(Ndecoup)]);

        coeffsErr = 1.96*std(coeffsRegLinMin.')/sqrt(Ndecoup);

        % save
        F0_err(ks) = coeffsErr(1);
        coeffAmpl_err(ks) = coeffsErr(2);
        coeffPos_err(ks) = coeffsErr(3);
        coeffTemp_err(ks) = coeffsErr(4);

        %% reg lin pos par heure

        coeffsPos = nan(1, 24);

        for kh = 1:24
            indexesHeure = (Tridge >= (kh-1)*3600) & (Tridge < kh*3600);
            coeffsRegLin0 = [ones(size(Fridge(indexesHeure))); Yridge(indexesHeure) - mean(Yridge(indexesHeure))].' \ Fridge(indexesHeure).';
            coeffsPos(kh) = coeffsRegLin0(2);
        end


        fprintf('\nReg. lin. pos :\ncoeff = %.2f +- %.2f Hz/m\n', [mean(coeffsPos), 1.96*std(coeffsPos)/sqrt(24)]);



        %% averaging

        if ~plotFA
            continue
        end

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

    end

    %%
    T = table(Ksimu, F0, F0_err, coeffAmpl, coeffAmpl_err, coeffPos, coeffPos_err, coeffTemp, coeffTemp_err);
    disp(' ');
    disp(ridgeFolder);
    disp(' ');
    disp(T);
%     save(fullfile(ridgeFolderPath, ridgeFolder, 'linreg.mat'), 'T');

end

