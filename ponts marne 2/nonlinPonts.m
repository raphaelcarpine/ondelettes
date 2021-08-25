clear all

removeNan = true;

plotTemp = 0;
plotFourier = 1;
plotRidge = 1;
plotHist = 1;
plotK = 1;
plotFA = 1;

projection = 5; % 1: mode1, 5: mode 5

%% data

% dataFolder = 'C:\Users\carpine\Documents\projets\ponts marne\reprise operations 2021\donnees'; % dossier où les fichier .csv sont
% dataFileName = 'esbly_1005.mat';
% dataFilePath = fullfile(dataFolder, dataFileName);

% dataFilePath = choixData();

bridges = 1:10;
bridges = 0;



Xtotal = {};
Ytotal = {};
stdYtotal = {};
Ktotal = {};
Namestotal = {};

if length(bridges) > 1
    [initWaitBar, updateWaitBar, closeWaitBar] = getWaitBar(length(bridges), 'windowTitle', 'Computing ridges',...
        'progressStringFcn', @(k) sprintf('(%u/%u)', [k, length(bridges)]), 'minimumBeepTime', 0);
else
    initWaitBar = @() 0;
    updateWaitBar = @(~, ~) 0;
    closeWaitBar = @() 0;
end
initWaitBar();

for kbridge = 1:length(bridges)
    Kbridge = bridges(kbridge);
    dataFilePath = choixData(Kbridge);
    
    load(dataFilePath);
    
    [~, dataFileName] = fileparts(dataFilePath);
    
    updateWaitBar(kbridge, strrep(dataFileName, '_', ' '));
    
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    disp(dataFileName);
    disp(startDate);
    % T°
    if plotTemp
        [TemperatureTime, TemperatureTemp] = getTemperature(startDate);
        figure('Name', dataFileName);
        plot(TemperatureTime, TemperatureTemp);
        ylabel('Temperature [°C]');
    end
    
    X = X.';
    T = T.';
    
    %% enlever donnees redondantes (double t à cause du mode transmit/log)
    
    [X, T] = removeRedundantData(X, T);
    
    
    %% enlever nan debut fin
    
    if removeNan
        [X, T] = removeNanSignal(X, T);
    end
    
    %% select channels
    % orthogonal projection
    
    if projection == 1
        ChanCenter1 = find(strcmp(channelNames, '40198:ch3'));
        ChanCenter2 = find(strcmp(channelNames, '29279:ch3'));
        X = X(ChanCenter1, :) + X(ChanCenter2, :);
    elseif projection == 5
        ChanSW = find(strcmp(channelNames, '40197:ch3'));
        ChanNW = find(strcmp(channelNames, '40199:ch3'));
        ChanSE = find(strcmp(channelNames, '29278:ch3'));
        ChanNE = find(strcmp(channelNames, '29280:ch3'));
        X = X(ChanSW, :) - X(ChanNW, :) - X(ChanSE, :) + X(ChanNE, :);
    else
        error();
    end
    
    
    %% CWT
    
    MotherWavelet = 'morlet';
    if projection == 1
        Q = sqrt(pi/0.01)/3;
        fmin = 2.1-1;
        fmax = 2.1+1;
%         fmin = 0.5;
%         fmax = 3.5;
    elseif projection == 5
        Q = sqrt(pi/0.01)/3;
        fmin = 4;
        fmax = 7;
%         fmin = 1; % test bruit
%         fmax = 3;
    end
    ct = 3;
    
    if 0
        fig = figure;
        ax = axes(fig);
        plts = plot(T, X);
        xlabel(ax, 'Time [s]');
        ylabel(ax, 'Acceleration [m/s²]');
        legend('40198:ch3 + 29279:ch3');
        
        
        NbMaxRidges = 10;
        NbMaxParallelRidges = 1;
        MaxSlopeRidge = inf;
        
        WaveletMenu('WaveletPlot', plts, 'MotherWavelet', MotherWavelet, 'Q', Q, 'fmin', fmin, 'fmax', fmax, 'RemoveMean', true,...
            'AutocorrelationMode', false, 'AutocorrelationMaxLag', 5/(2*pi*2.2*0.01), 'FourierScale', 'log',...
            'MaxRidges', NbMaxRidges, 'MaxParallelRidges', NbMaxParallelRidges, 'MaxSlopeRidge', MaxSlopeRidge);
        return
    end
    
    %%
    
    X = X - mean(X);
    
    % fourier
    [f, Fx] = fourierTransform(T, X, 'Averaging', true, 'AveragingNb', 100);
    if plotFourier
        figure('Name', dataFileName);
        plot(f, Fx);
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [m/s²]');
        xlim([0, 20]);
    end
    
    [~, kfMax] = max(Fx);
    fprintf('Mode freq: %.2f Hz\n\n', f(kfMax));
    
    
    % cwt
    freqs = linspace(fmin, fmax, 300);
    
    CWT = WvltComp(T, X, freqs, Q, 'MotherWavelet', MotherWavelet);
    [~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
    
    Tridge = T(T >= T(1)+ct*DeltaT(fmin) & T <= T(end)-ct*DeltaT(fmin));
    CWTridge = [nan(1, length(Tridge)); CWT(:, T >= T(1)+ct*DeltaT(fmin) & T <= T(end)-ct*DeltaT(fmin)); nan(1, length(Tridge))];
    freqs = [nan, freqs, nan];
    [~, Fridge] = max(abs(CWTridge), [], 1);
    [Fridge, Aridge] = localMax3Points(freqs([Fridge-1; Fridge; Fridge+1]),...
        CWTridge([Fridge-1; Fridge; Fridge+1] + [1;1;1] * (0:size(CWTridge, 2)-1)*size(CWTridge, 1)));
    Aridge = abs(Aridge);
    
    
    %%
    
    if 0
        figure('Name', dataFileName);
        plot(Aridge, Fridge);
        xlabel('Amplitude [m/s²]');
        ylabel('Frequency [Hz]');
        
        ScatterAveragingMenu;
        return
    end
    
    if plotRidge
        figure('Name', dataFileName);
        plot(Aridge, Fridge);
        set(gca, 'XScale', 'log');
        xlabel('Amplitude [m/s²]');
        ylabel('Frequency [Hz]');
    end
    
    Aridge = Aridge(~isnan(Fridge));
    Tridge = Tridge(~isnan(Fridge));
    Fridge = Fridge(~isnan(Fridge));
    
    if plotHist
        figure('Name', dataFileName);
        hist(log10(Aridge), 100);
        set(gca,'Xtick',-5:1);
        set(gca,'Xticklabel',10.^get(gca,'Xtick'));
        xlabel('Amplitude [m/s²]');
    end
    
    [X, Y, stdY, K, Xlims] = averagingScatter3(Aridge, Fridge, 0.1, 'log', Tridge, 100);
    
    Kmin = 3;
    Y(K < Kmin) = nan;
    stdY(K < Kmin) = nan;
    
    if plotK
        figure('Name', dataFileName);
        bar(log10(X), K);
        set(gca,'Xtick',-5:1);
        set(gca,'Xticklabel',10.^get(gca,'Xtick'));
        % set(gca, 'XScale', 'log');
        xlabel('Amplitude [m/s²]');
    end
    
    if plotFA
        figure('Name', dataFileName);
        err = 1.96*stdY./sqrt(K);
        errorbar(X, Y, err);
        set(gca, 'XScale', 'log');
        xlabel('Amplitude [m/s²]');
        ylabel('Frequency [Hz]');
    end
    
    %%
    Xtotal{end+1} = X;
    Ytotal{end+1} = Y;
    stdYtotal{end+1} = stdY;
    Ktotal{end+1} = K;
    Namestotal{end+1} = strrep(dataFileName, '_', ' ');

end

closeWaitBar();

if length(Xtotal) > 1
    figure;
    for k = 1:length(Xtotal)
        plot(Xtotal{k}, Ytotal{k});
        hold on
    end
    set(gca, 'XScale', 'log');
    legend(Namestotal);
end

selectLine(gca);

