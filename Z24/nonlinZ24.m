clear all
% close all

filtering = false;
testCWT = false;
averageScale = 'lin';
plotDiscardedFA = true;

pltFFT = 0;
plotRidge = 0;
plotHist = 0;
plotK = 0;
plotFA = 1;

vibration = 'forced';
measures = 1:17; % 1-17, or 1:17

%% loop

plots = [pltFFT, plotRidge, plotHist, plotK, plotFA];
if length(measures) > 1 && sum(plots) > 1
    fprintf('%u plots for each measure, continue? ', sum(plots));
    cont = input('plot signal? ', 's');
    if ~ismember({cont}, {'', 'y', 'yes'})
        return
    end
end

for measure  = measures
    
    fprintf('measure %02u\n', measure);
    
    %% mode 1
    
    config = 5; % only config for mode 1
    [X, t, labels] = getDataZ24(measure, config, vibration);
    
    I = false(size(labels));
    for ki = 1:length(I)
        if labels{ki}(1) ~= 'R' && labels{ki}(1) ~= 'D' && labels{ki}(end) == 'V'
            I(ki) = true;
        end
    end
    
    x = sum(X(I, :), 1);
    x = x - mean(x);
    
    % check
    if sum(I) ~= 15
        warning(sprintf('Nb ch. = %u =/= 15', sum(I)));
    end
    
    %% filtering (high pass)
    
    F0 = 3;
    Nfiltering = 5;
    
    if filtering
        x = butterworthFilter(t, x, F0, 'high', Nfiltering);
    end
    
    %% fft
    
    [f, Fx] = fourierTransform(t, x, 'Averaging', true, 'AveragingNb', 20);
    
    if pltFFT
        figure('Name', sprintf('measure %02u, config %u, fft mode 1', [measure, config]));
        fftPlt = plot(f, Fx);
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [m/s²]');
        xlim([0, 10]);
        
        % peak picking
        [fmax, zeta] = PeakPickingMax('Line', fftPlt, 'QuadraticMax', false,...
            'QuadraticFFT', false, 'PlotPeak', true);
    else
        [fmax, zeta] = PeakPickingMax('Freqs', f, 'FFT', Fx, 'QuadraticMax', false,...
            'QuadraticFFT', false, 'PlotPeak', false);
    end
    fprintf('  Fmax = %.2f Hz ; zeta = %.2f %%\n', [fmax, 100*zeta]);
    
    %% CWT
    
    % CWT
    ct = 3;
    cA = 3;
    fmin = fmax - 1;
    fmax = fmax + 1;
    Q = sqrt(pi/zeta)/cA;
    MotherWavelet = 'morlet';
    
    if testCWT % wavelet menu
        fig = figure('Name', sprintf('measure %02u, config %u', [measure, config]));
        plt2 = plot(t, x);
        xlabel('Time [s]');
        ylabel('Acceleration [m/s²]');
        selectLine(gca);
        
        WaveletMenu('WaveletPlot', plt2, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
            'MotherWavelet', MotherWavelet, 'MaxParallelRidges', 1, 'MaxSlopeRidge', inf);
        continue
    end
    
    % one ridge extraction
    freqs = linspace(fmin, fmax, 300);
    
    CWT = WvltComp(t, x, freqs, Q, 'MotherWavelet', MotherWavelet);
    [~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
    
    Tridge = t(t >= t(1)+ct*DeltaT(fmin) & t <= t(end)-ct*DeltaT(fmin));
    CWTridge = [nan(1, length(Tridge)); CWT(:, t >= t(1)+ct*DeltaT(fmin) & t <= t(end)-ct*DeltaT(fmin)); nan(1, length(Tridge))];
    freqs = [nan, freqs, nan];
    [~, Fridge] = max(abs(CWTridge), [], 1);
    [Fridge, Aridge] = localMax3Points(freqs([Fridge-1; Fridge; Fridge+1]),...
        CWTridge([Fridge-1; Fridge; Fridge+1] + [1;1;1] * (0:size(CWTridge, 2)-1)*size(CWTridge, 1)));
    Aridge = abs(Aridge);
    
    Aridge = Aridge(~isnan(Fridge));
    Tridge = Tridge(~isnan(Fridge));
    Fridge = Fridge(~isnan(Fridge));
    
    % one ridge plot
    if plotRidge
        figure('Name', sprintf('ridge, Q = %.1f (measure %02u, config %u)', [Q, measure, config]));
        plot(Aridge, Fridge);
        xlabel('Amplitude [m/s²]');
        ylabel('Frequency [Hz]');
    end
    
    % hist amplitudes
    if plotHist
        figure('Name', sprintf('ridge hist, Q = %.1f (measure %02u, config %u)', [Q, measure, config]));
        hist(log10(Aridge), 30);
        set(gca,'Xtick',-5:1);
        set(gca,'Xticklabel',10.^get(gca,'Xtick'));
        xlabel('Amplitude [m/s²]');
    end
    
    % averaging
    if strcmp(averageScale, 'lin')
        incrementA = 0.001;
    else
        incrementA = 0.1;
    end
    [X, Y, stdY, K, Xlims] = averagingScatter3(Aridge, Fridge, incrementA, averageScale, Tridge, 100);
    
    Kmin = 3;
    Y0 = Y;
    stdY0 = stdY;
    Y(K < Kmin) = nan;
    stdY(K < Kmin) = nan;
    
    if plotK
        figure('Name', sprintf('K averaging, Q = %.1f (measure %02u, config %u)', [Q, measure, config]));
        if strcmp(averageScale, 'lin')
            bar(X, K);
        else
            bar(log10(X), K);
            set(gca,'Xtick',-5:1);
            set(gca,'Xticklabel',10.^get(gca,'Xtick'));
        end
        % set(gca, 'XScale', 'log');
        xlabel('Amplitude [m/s²]');
    end
    
    if plotFA
        figure('Name', sprintf('ridge average (with 95%% confidence interval), Q = %.1f (measure %02u, config %u)', [Q, measure, config]));
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
        xlabel('Amplitude [m/s²]');
        ylabel('Frequency [Hz]');
        if length(measures) > 1
            if strcmp(averageScale, 'lin')
                xlim([0 0.045]);
            else
                xlim([1e-5 0.045]);
            end
        end
        ylim([3 4.5]);
    end
    
    drawnow;
    
    
end




