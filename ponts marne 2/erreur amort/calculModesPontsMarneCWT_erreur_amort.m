clear all

savePath = 'ponts marne 2\erreur amort\save';

MotherWavelet = 'cauchy';

% mode
manualMode = 0;
saveResults = 1;
saveResults = saveResults & ~manualMode;

% affichage
pauseModes = 0; % pause entre les modes
plotCWT = 0;
plotRidgeExtract = 0;
plotShapes = 0;
plotShapesTime = 0;
plotTemp = 0;

% filtrage passe haut
filtrage = 1;
fc_filtre = 3.5; % freq coupure
fmin_filtrage = 5; % min freq propre avec filtre

% séparation signal
Nsep = 30;

% pont et mode
bridge = 6;
Kf = 0;

% Meff
Meff0 = 1;

%% data

dataFilePath = choixData(bridge);

load(dataFilePath);

[~, dataFileName] = fileparts(dataFilePath);
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

% mise en forme
X = X.';
T = T.';

[X, T] = removeRedundantData(X, T);

[X, T] = removeNanSignal(X, T);

X = -X; % capteurs vers le bas
X = X - mean(X, 2); % moyenne


% tri channels par ordre croissant
[channelNames, I] = sort(channelNames);
X = X(I, :);
I = 1:length(I);
for kc = 1:length(channelNames) % ch2 à la fin
    if strcmp(channelNames{kc}(end-2:end), 'ch2')
        channelNames = [channelNames(1:kc-1), channelNames(kc+1:end), channelNames(kc)];
        I = [I(1:kc-1), I(kc+1:end), I(kc)];
        break
    end
end
X = X(I, :);

Nt = length(T);
dt = (T(end) - T(1))/(Nt-1);
Ttot = Nt*dt;

% découpage sous signaux
Ntsep = floor(Nt/Nsep); % longueur sous-intervalles
Isep = cell(1, Nsep); % sous-intervalles
for ksep = 1:Nsep
    Isep{ksep} = (ksep-1)*Ntsep + (1:Ntsep);
end
Tsep = Ntsep*dt;

% deformees modales
dimensionsShapes2;

%% calcul cwt

% freqs approx CWT
if bridge == 6
    Freqs0 =      [2.08 2.18 2.84 2.94 5.47 6.52 7.85 15.04 16.74 21.13];
    Damps0 = 0.01*[0.9  0.9  0.8  1.0  1.8  3.7  5.2  2.4   1.0   1.0  ];
    Freqs0bis = [Freqs0, 19.31];
else
    error(' ');
end

if isempty(Kf) || (length(Kf) == 1 && Kf == 0)
    Kf = 1:length(Freqs0);
end

% save
Freqs = nan(length(Kf), Nsep);
Damps = nan(length(Kf), Nsep);
QualityFactors = nan(length(Kf), Nsep);
Tinit = nan(length(Kf), Nsep);
Shapes = nan(size(X, 1), length(Kf), Nsep);
PbCalculRidge = false(length(Kf), Nsep);

% waitbar
[initWaitBar, updateWaitBar, closeWaitBar] =...
    getWaitBar(length(Kf)*Nsep, 'windowTitle', 'Extraction modes', 'displayTime', 0);
initWaitBar();

for ikf = 1:length(Kf)
    kf = Kf(ikf);
    
    freq0 = Freqs0(kf);
    damp0 = Damps0(kf);
    
    % pole
    poleMode = 2i*pi*freq0 - 2*pi*damp0*freq0/sqrt(1-damp0^2);
    alphap = sqrt((-4/real(poleMode) + 2/abs(imag(poleMode)))/Tsep);
    Topt = -1.4326/real(poleMode);
    
    Df = sort(abs(Freqs0bis - freq0));
    Df = Df(2);
    ct = 3;
    cf = 5;
    Qmin =  getBoundsQ2(freq0, Df, 0, [0 0], [0 0], ct, ct, cf, MotherWavelet);
    
    % gestion cas modes proches
    Nsv = 1;
    
    % filtrage
    if filtrage && freq0 > fmin_filtrage
        Xfiltr = butterworthFilter(T, X, fc_filtre, 'high', 10);
        fprintf('filtrage (fc = %.2fHz)\n', fc_filtre);
    else
        Xfiltr = X;
    end
    
    fmin = freq0-Df/2;
    fmax = freq0+Df/2;
    NbFreqs = 100;
    Q = Qmin;
    [~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
    Dt_edge = ct*DeltaT(freq0);
    
    % display
    fprintf('\n~~~~~~ mode %d (%d/%d) ~~~~~~\n', [kf, ikf length(Kf)]);
    fprintf('f = %.2fHz, z = %.2f%%, 1/mu = %.2fs\n', [freq0, 100*damp0, -1/real(poleMode)]);
    fprintf('Q = %.1f, tri = %.2fs, Topt = %.2fs, Tsep = %.1fs\n', [Q, Dt_edge, Topt, Tsep]);
    fprintf('alpha = %.2f (= %.2f t=tri, = %.2f t=tri+Topt)\n', [alphap, alphap*exp(-real(poleMode)*Dt_edge), alphap*exp(-real(poleMode)*(Dt_edge+Topt))]);
    fprintf('cov z = %.2g\n', Meff0 * exp(-real(poleMode*Dt_edge)) * alphap);
    
    for ksep = 1:Nsep
        
        updateWaitBar(nan, sprintf('mode %d (%d/%d)', [kf, ikf, length(Kf)]));
        
        if manualMode
            fig = figure;
            plt = plot(T(Isep{ksep}), Xfiltr(:, Isep{ksep}));
            WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'XLimRidge', [0, Dt_edge + Topt],...
                'AutocorrelationMode', true, 'AutocorrelationSVDMode', true, 'AutocorrelationFourierSVDMode', true,...
                'AutocorrelationMaxLag', 3*Dt_edge + Topt, 'AutocorrelationNsvd', Nsv, 'FourierScale', 'lin',...
                'MaxRidges', 1, 'MaxParallelRidges', inf, 'RealShapePlot', shapePlotBridge,...
                'AnimatedShapePlot', shapePlotBridgeAnim, 'MotherWavelet', MotherWavelet);
            waitfor(fig);
            continue
        end
        
        % calcul cross-corr
        Rx = crossCorrelation(Xfiltr(:, Isep{ksep}), round((2.5*Dt_edge + Topt)/dt), 'unbiased');
        tRx = dt * (0:size(Rx, 3)-1);
        
        % calcul CWT
        [SVrx, SVvectrx] = svdCWT(tRx, Rx, linspace(fmin, fmax, NbFreqs), Q, Nsv,...
            'MotherWavelet', MotherWavelet);
        
        % calcul ridge
        [t_r, freqs_r, shapes_r, amplitudes_r] = getModesCrossCorr(tRx, SVrx, SVvectrx,...
            Q, fmin, fmax, NbFreqs, Nsv, 'NbMaxRidges', 1, 'NbMaxParallelRidges', inf,...
            'XLimRidge', [0, Dt_edge + Topt], 'MotherWavelet', MotherWavelet, 'StopWhenIncreasing', false);
        
        % cas pas de ridge
        if isempty(t_r{Nsv})
            warning('problème extraction ridge (pas de ridge)');
            PbCalculRidge(ikf, ksep) = true;
            continue
        end
        
        % conversion
        t_r = t_r{Nsv}{1};
        freqs_r = freqs_r{Nsv}{1};
        shapes_r = shapes_r{Nsv}{1};
        amplitudes_r = amplitudes_r{Nsv}{1};
        
        % cas d'aret du ridge avant la limite
        if abs(t_r(end) - (Dt_edge + Topt)) > 3*dt
            fprintf('(arrêt ridge avant borne) ');
            PbCalculRidge(ikf, ksep) = true;
        end
        % cas début ridge après début
        if abs(t_r(1) - ct*DeltaT(freqs_r(1))) > 3*dt
            fprintf('(début ridge après borne) ');
            PbCalculRidge(ikf, ksep) = true;
        end
        
        % calcul moyennes
        shape_moy = mean(shapes_r, 2);
        shape_moy = shape_moy * sign(real(shape_moy(8))); %orientation
        freq_moy = mean(freqs_r);
        coeffsRegLin = log(amplitudes_r) / [ones(size(t_r)); t_r];
        lambda_reg = coeffsRegLin(2);
        
        % parametres modaux
        poleMode2 = 2i*pi*freq_moy + lambda_reg;
        freq2 = abs(poleMode2) / (2*pi);
        damp2 = -real(poleMode2) / abs(poleMode2);
        
        % print
        fprintf('%.2f, ', 100*damp2);
        
        % save
        Freqs(ikf, ksep) = freq2;
        Damps(ikf, ksep) = damp2;
        Shapes(:, ikf, ksep) = shape_moy;
        QualityFactors(ikf, ksep) = Q;
        Tinit(ikf, ksep) = Dt_edge;
        
        % waitbar
        updateWaitBar();
        
        % plot
        if plotRidgeExtract
            % freq
            fig = figure;
            fig.Position(1) = fig.Position(1) - fig.Position(3)/2;
            plot(t_r, freqs_r);
            yline(freq_moy, '--r');
            xline(Dt_edge + Topt);
            set(gca, 'ylim', freq_moy*(1 + 0.1*[-1 1]));
            % damp
%             Rxx = nan(size(Rx, 1), size(Rx, 3));
%             for kdof = 1:size(Rx, 1)
%                 Rxx(kdof, :) = Rx(kdof, kdof, :);
%             end
            fig = figure;
            fig.Position(1) = fig.Position(1) + fig.Position(3)/2;
%             plot(tRx, Rxx);
%             hold on
            plot(t_r, amplitudes_r, 'LineWidth', 2);
            hold on
            plot(t_r, exp(coeffsRegLin * [ones(size(t_r)); t_r]), '--r', 'LineWidth', 2);
            xline(Dt_edge + Topt);
            set(gca, 'yscale', 'log');
            set(gca, 'ygrid', 'on');
        end
        
        if plotCWT
            plt = WvltPlot2(tRx, linspace(fmin, fmax, NbFreqs), SVrx{1}, 'abs', Q, ct, MotherWavelet, 'log');
            fig = gcf;
            fig.Position(1) = fig.Position(1) - fig.Position(3)/2;
        end
        
        if plotShapes
            figName = sprintf('f = %.2fHz, z = %.2f%%, I = %.2f %%', [freq0, 100*damp0, 100*nonPropIndex(shape_moy)]);
            fig = shapePlotBridge(real(shape_moy), figName);
            fig.Position(1) = fig.Position(1) - fig.Position(3)/2;
            fig.Position(2) = fig.Position(2) - fig.Position(4) - 10;
            
            %         fig = shapePlotBridgeAnim(shape_moy, figName);
            %         fig.Position(1) = fig.Position(1) + fig.Position(3)/2;
            %         fig.Position(2) = fig.Position(2) - fig.Position(4) - 10;
        end
        
        if plotShapesTime
            fig = figure;
            fig.Position(1) = fig.Position(1) - fig.Position(3)/2;
            plot(t_r, real(shapes_r));
            xlabel('t');
            ylabel('Re');
            
            fig = figure;
            fig.Position(1) = fig.Position(1) + fig.Position(3)/2;
            plot(t_r, imag(shapes_r));
            xlabel('t');
            ylabel('Im');
        end
        
        if pauseModes
            input('continue?');
        end
        
        close all
        
    end
    
    % results
    fprintf('\n\nfmoy = %.2fHz, std = %.2fHz\n', [mean(Freqs(ikf, :)), std(Freqs(ikf, :))]);
    fprintf('zmoy = %.2f%%, std = %.2f%%\n', 100*[mean(Damps(ikf, :)), std(Damps(ikf, :))]);
    Meff = std(Damps(ikf, :))/damp0 * exp(real(poleMode)*Dt_edge) / alphap;
    fprintf('Meff = %.2f\n', Meff);
end
closeWaitBar();


if saveResults
    saveFile = fullfile(savePath, ['modes_', dataFileName, '_Nsep', num2str(Nsep)]);
    save(saveFile, 'Freqs', 'Damps', 'QualityFactors', 'Tinit', 'Shapes', 'PbCalculRidge', 'Tsep');
    disp('saved');
end




















