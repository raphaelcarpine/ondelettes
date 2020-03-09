load('mur silvia\modesFourier\ModesLMS.mat');


%%
% affichage et methode de calcul
singleRidgeMode = true;

verb = true;
plotTemporel = true;
saveFiles = true;
progressBar = true;

% choix precision
ct = 3;
cf = 5;

if verb
    progressBar = false;
end

%% moyennes et ecarts types

% nombre de modes pour P0, P6 et P7
nbModes = [3, 4, 3];

% variables de sauvegarde

for indp = 1:3
    for mode = 1:nbModes(indp)
        allFreqs{indp}{mode} = [];
        allShapes{indp}{mode} = [];
        allDamps{indp}{mode} = [];
    end
end


% meanFreqs = {[], [], []};
% meanShapes = {[], [], []};
% meanDamps = {[], [], []};
% nbTransients = {[], [], []};
% stdFreqs = {[], [], []};
% stdShapes = {[], [], []};
% stdDamps = {[], [], []};

%% data

P = [0, 6, 7];

% données de l'algo polymax
Freqs = {[], [], []};
Damps = {[], [], []};
for p = 1:3
    mode = 1;
    while mode <= size(ModesLMS, 2) && ~isempty(ModesLMS(p, mode).freq)
        Freqs{p}(end+1) = ModesLMS(p, mode).freq;
        Damps{p}(end+1) = ModesLMS(p, mode).damping;
        mode = mode+1;
    end
end


% ModesTransients = {{[1; 26], [], [1, 2, 3; 26, 3, 3]},... % P0
%     {[1, 2, 3; 3, 1.8, 4.3], [1, 2, 3; 1.8, 1.4, 1.4]},... % P6
%     {[1, 2; 17.4, 2.4], [2; 2.5], [1; 17.4]}}; % P7

TransientsModes = {[],... % P0
    [],... % P6
    []}; % P7

TransientsNumbers = {[],... % P0
    [],... % P6
    []}; % P7

TransientsTimes = {[],... % P0
    [],... % P6
    []}; % P7

TransientsDeltaF = {[],... % P0
    [],... % P6
    []}; % P7


%P0T0
TransientsModes{1} = [TransientsModes{1}, 1];
TransientsNumbers{1} = [TransientsNumbers{1}, 0];
TransientsTimes{1} = [TransientsTimes{1}, [700.2; 701.5]];
TransientsDeltaF{1} = [TransientsDeltaF{1}, 2.7];

%P0T1
TransientsModes{1} = [TransientsModes{1}, 2, 3];
TransientsNumbers{1} = [TransientsNumbers{1}, 1, 1];
TransientsTimes{1} = [TransientsTimes{1}, [1218.1, 1218.1; 1223, 1219.3]];
TransientsDeltaF{1} = [TransientsDeltaF{1}, 2.8, 2.8];

%P0T2
TransientsModes{1} = [TransientsModes{1}, 2, 3];
TransientsNumbers{1} = [TransientsNumbers{1}, 2, 2];
TransientsTimes{1} = [TransientsTimes{1}, [1267.35, 1267.35; 1269.95, 1269.95]];
TransientsDeltaF{1} = [TransientsDeltaF{1}, 2.8, 2.8];

%P0T3
TransientsModes{1} = [TransientsModes{1}, 1, 2];
TransientsNumbers{1} = [TransientsNumbers{1}, 3, 3];
TransientsTimes{1} = [TransientsTimes{1}, [1450, 1450; 1452.5, 1453]];
TransientsDeltaF{1} = [TransientsDeltaF{1}, 2.7, 2.8];

%P6T1
TransientsModes{2} = [TransientsModes{2}, 1, 2, 3, 4]; % 3 à vérifier
TransientsNumbers{2} = [TransientsNumbers{2}, 1, 1, 1, 1];
TransientsTimes{2} = [TransientsTimes{2}, [238, 238, 239, 238; 243, 241, 243, 239.5]]; % t0 m3
TransientsDeltaF{2} = [TransientsDeltaF{2}, 4.7, 1.5, 1.5, 4.9];

%P6T2
TransientsModes{2} = [TransientsModes{2}, 1];
TransientsNumbers{2} = [TransientsNumbers{2}, 2];
TransientsTimes{2} = [TransientsTimes{2}, [278.6; 282.6]];
TransientsDeltaF{2} = [TransientsDeltaF{2}, 4.7];

%P6T3
TransientsModes{2} = [TransientsModes{2}, 3];
TransientsNumbers{2} = [TransientsNumbers{2}, 3];
TransientsTimes{2} = [TransientsTimes{2}, [444; 447]];
TransientsDeltaF{2} = [TransientsDeltaF{2}, 1.4];

%P7T1
TransientsModes{3} = [TransientsModes{3}, 3];
TransientsNumbers{3} = [TransientsNumbers{3}, 1];
TransientsTimes{3} = [TransientsTimes{3}, [6.72; 7.83]];
TransientsDeltaF{3} = [TransientsDeltaF{3}, 3.4];

%P7T2
TransientsModes{3} = [TransientsModes{3}, 2, 3];
TransientsNumbers{3} = [TransientsNumbers{3}, 2, 2];
TransientsTimes{3} = [TransientsTimes{3}, [36.81, 36.81; 38.5, 39]];
TransientsDeltaF{3} = [TransientsDeltaF{3}, 2.3, 3.4];

%P7T3
TransientsModes{3} = [TransientsModes{3}, 2, 3];
TransientsNumbers{3} = [TransientsNumbers{3}, 3, 3];
TransientsTimes{3} = [TransientsTimes{3}, [61.7, 61.7; 63.2, 63.2]];
TransientsDeltaF{3} = [TransientsDeltaF{3}, 2.3, 3.4];

%P7T4
TransientsModes{3} = [TransientsModes{3}, 1];
TransientsNumbers{3} = [TransientsNumbers{3}, 4];
TransientsTimes{3} = [TransientsTimes{3}, [267.5; 271.4]];
TransientsDeltaF{3} = [TransientsDeltaF{3}, 2.3];

%P7T5
TransientsModes{3} = [TransientsModes{3}, 3];
TransientsNumbers{3} = [TransientsNumbers{3}, 5];
TransientsTimes{3} = [TransientsTimes{3}, [285.8; 286.7]];
TransientsDeltaF{3} = [TransientsDeltaF{3}, 3.4];

%P7T6
TransientsModes{3} = [TransientsModes{3}, 1, 2];
TransientsNumbers{3} = [TransientsNumbers{3}, 6, 6];
TransientsTimes{3} = [TransientsTimes{3}, [286.7, 286.7; 289, 289]];
TransientsDeltaF{3} = [TransientsDeltaF{3}, 2.3, 2.3];

%P7T7
TransientsModes{3} = [TransientsModes{3}, 2, 3];
TransientsNumbers{3} = [TransientsNumbers{3}, 7, 7];
TransientsTimes{3} = [TransientsTimes{3}, [305.9, 305.9; 306.8, 307.5]];
TransientsDeltaF{3} = [TransientsDeltaF{3}, 6.2, 6.2];

%P7T8
TransientsModes{3} = [TransientsModes{3}, 1, 3];
TransientsNumbers{3} = [TransientsNumbers{3}, 8, 8];
TransientsTimes{3} = [TransientsTimes{3}, [324.7, 324.7; 326.5, 326]];
TransientsDeltaF{3} = [TransientsDeltaF{3}, 2.3, 3.4];

%P7T9
TransientsModes{3} = [TransientsModes{3}, 3];
TransientsNumbers{3} = [TransientsNumbers{3}, 9];
TransientsTimes{3} = [TransientsTimes{3}, [335.1; 336.3]];
TransientsDeltaF{3} = [TransientsDeltaF{3}, 3.4];

%P7T10
TransientsModes{3} = [TransientsModes{3}, 2, 3];
TransientsNumbers{3} = [TransientsNumbers{3}, 10, 10];
TransientsTimes{3} = [TransientsTimes{3}, [458.1, 458.1; 459.8, 459.8]];
TransientsDeltaF{3} = [TransientsDeltaF{3}, 6.2, 6.2];






%%

% barre progression
if progressBar
    Nprogress = length([TransientsModes{:}]);
    Kprogress = 0;
    progressBarFig = waitbar( Kprogress/Nprogress, '');
end


for indp = 2:3
    p = P(indp);
    freqs = Freqs{indp};
    damps = Damps{indp};
    
    TransientsModesP = TransientsModes{indp};
    TransientsNumbersP = TransientsNumbers{indp};
    TransientsTimesP = TransientsTimes{indp};
    TransientsDeltaFP = TransientsDeltaF{indp};
    
    if verb
        disp(' ');
        disp(['~~~~~~ P', num2str(p), ' ~~~~~~']);
        disp(' ');
    end
    
    for kridge = 1:length(TransientsNumbersP)
        transient = TransientsNumbersP(kridge);
        mode = TransientsModesP(kridge);
        
        if verb
            disp(['~~~ P', num2str(p), 'T', num2str(transient)]);
            disp(['~ mode', num2str(mode)]);
        end
        
        f = freqs(mode);
        Df = TransientsDeltaFP(kridge);
        damp = damps(mode);
        Dt = 1 / (damp * 2*pi*f);
        
        % chargement des donnees
        [t, X] = getData(p, 0);
        
        boundsT = TransientsTimesP(:, kridge);
        X = X(:, t>=boundsT(1) & t<boundsT(2));
        t = t(t>=boundsT(1) & t<boundsT(2));
        
        T = t(end) - t(1);
        
        % choix Q
        [Qmin, Qmax, Qz] = getBoundsQ(f, Df, Dt, T, ct, cf);
        if Qmin > min(Qmax, Qz)
            warning('Qmin > min(Qmax, Qz)');
        end
        
        Q = (Qmin + min(Qmax, Qz)) / 2;
        Q = Qmin;
        if verb
            disp(['Qmin = ', num2str(Qmin), ' ; Qmax = ', num2str(Qmax), ' ; Qz = ', num2str(Qz)]);
            disp(['Q : ', num2str(Q)]);
        end
        
        % calcul modes
        fmin = f - min(f/2, Df/2);
        fmax = f + min(f/2, Df/2);
        NbFreq = 300;
        
        MaxRidges = 1;
        MaxParallelRidges = 1;
        
        
        if singleRidgeMode
            [time, freq, shape, amplitude] = getModesSingleRidge(t, X, Q, fmin, fmax, NbFreq,...
                'NbMaxRidges', MaxRidges, 'NbMaxParallelRidges', MaxParallelRidges,...
                'ctLeft', ct, 'ctRight', ct);
        else
            ridges = {};
            for k = 1:9
                ridges{end+1} = RidgeExtract(t, X(k,:), Q, fmin, fmax, NbFreq,...
                    'NbMaxParallelRidges', MaxRidges, 'NbMaxRidges', MaxParallelRidges,...
                    'ctLeft', ct, 'ctRight', ct);
            end
            
            [time, freq, ~, shape, amplitude, errors, ridgesNumber] = getModes(ridges, 1);
        end
        
        if isempty(time)
            warning(['no ridge, P', num2str(p), 'T', num2str(transient), 'm', num2str(mode)]);
            if verb
                % test
                [ttest, Xtest] = getData(p, 0);
                boundsT = TransientsTimesP(:, kridge);
                figure;
                plts = plot(ttest, Xtest(1:9,:));
                plts = transpose(plts);
                
                WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
                    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges, 'CtEdgeEffects', ct,...
                    'XLim', boundsT);
                
                input('continue ?');
            end
            continue
        end
        
        % threshold noise
        thresholdNoise = nan(1, 9);
        sumWvltF2 = 0;
        for k = 1:9
            WvltF = WvltComp(t, X(k, :), f, Q);
            thresholdNoise(k ) = exp( mean( log( abs( WvltF))));
            sumWvltF2 = sumWvltF2 + WvltF.^2;
        end
        thresholdNoiseTotal = exp( mean( 1/2 * log( abs( sumWvltF2))));
        
        % threshold
        timeThreshold = time{1};
        freqThreshold = freq{1};
        shapeThreshold = shape{1};
        amplitudeThreshold = amplitude{1};
        if singleRidgeMode
            timeThreshold = timeThreshold( abs(amplitudeThreshold) >= thresholdNoiseTotal);
            freqThreshold = freqThreshold( abs(amplitudeThreshold) >= thresholdNoiseTotal);
            shapeThreshold = shapeThreshold(:, abs(amplitudeThreshold) >= thresholdNoiseTotal);
            amplitudeThreshold = amplitudeThreshold( abs(amplitudeThreshold) >= thresholdNoiseTotal);
        else
            % à faire
        end
        
        % moyenne
        meanFreq = mean(freqThreshold);
        shapeNoNan = shapeThreshold;
        shapeNoNan(isnan(shapeNoNan)) = 0;
        meanShape = mean(shapeNoNan, 2);
        
        funReg = @(param) log(param(1)) - param(2)*(timeThreshold-timeThreshold(1)) - log(abs(amplitudeThreshold));
        options = optimoptions(@lsqnonlin, 'Display', 'off');
        A0lambda = lsqnonlin(funReg, [abs(amplitudeThreshold(1)), damp], [], [],...
            options);
        A0 = A0lambda(1);
        lambda = A0lambda(2);
        zeta = lambda / (2*pi*meanFreq);
        
        % sauvegarde
        allFreqs{indp}{mode}(end+1) = meanFreq;
        allDamps{indp}{mode}(end+1) = zeta;
        allShapes{indp}{mode} = [allShapes{indp}{mode}; transpose(meanShape)];
        
        % comparaison fourier
        errorFreq = abs(meanFreq - ModesLMS(indp, mode).freq) / ModesLMS(indp, mode).freq;
        
        shapeF = ModesLMS(indp, mode).shape;
        
        errorShape = meanShape - shapeF;
        errorShape = sqrt( errorShape'*errorShape / (shapeF'*shapeF));
        
        % modal assurance criterion
        mac = abs(meanShape'*shapeF)^2 / ((meanShape'*meanShape) * (shapeF'*shapeF));
        
        if verb
            disp(['freq : ', num2str(meanFreq), ' ; freq lms : ', num2str(ModesLMS(indp, mode).freq), ' ; error : ', num2str(100*errorFreq), '%']);
            disp(['MAC : ', num2str(100*mac), '%, ', 'shape error : ', num2str(100*errorShape), '%']);
            disp(['amort : ', num2str(100*zeta), '% ; amort lms : ', num2str(100*ModesLMS(indp, mode).damping), '%']);
            disp(['I : ', num2str(100*nonPropIndex(meanShape)), '%',...
                ' ; I lms : ', num2str(100*nonPropIndex(ModesLMS(indp, mode).shape)), '%']);
        end
        
        % plots temporels
        figs = [];
        if plotTemporel && verb
            figure;
            hold on
            plot(time{1}, (angle(shape{1}*exp(-1i*pi/2)) + pi/2) * 180/pi, ':');
            ax = gca;
            ax.ColorOrderIndex = 1;
            plot(timeThreshold, (angle(shapeThreshold*exp(-1i*pi/2)) + pi/2) * 180/pi);
            plot(timeThreshold, zeros(size(timeThreshold)), 'black--');
            plot(timeThreshold, 180*ones(size(timeThreshold)), 'black--');
            ylim([-90, 270]);
            xlabel('t');
            ylabel('arg(T)');
            
            figure;
            hold on
            plot(time{1}, real(shape{1}), ':');
            ax = gca;
            ax.ColorOrderIndex = 1;
            plot(timeThreshold, real(shapeThreshold));
            xlabel('t');
            ylabel('Re(T)');
            figure;
            hold on
            plot(time{1}, imag(shape{1}), ':');
            ax = gca;
            ax.ColorOrderIndex = 1;
            plot(timeThreshold, imag(shapeThreshold));
            xlabel('t');
            ylabel('Im(T)');
            
            figure;
            hold on
            plot(time{1}, freq{1}, ':');
            ax = gca;
            ax.ColorOrderIndex = 1;
            plot(timeThreshold, freqThreshold);
            xlabel('t');
            ylabel('f');
            
            fig = figure;
            hold on
            plot(time{1}, abs(amplitude{1}), ':');
            ax = gca;
            ax.ColorOrderIndex = 1;
            plot(timeThreshold, abs(amplitudeThreshold));
            plot(time{1}, A0*exp(-lambda*(time{1}-time{1}(1))), 'r--');
            plot(time{1}, thresholdNoiseTotal * ones( size( time{1})));
            set(ax, 'YScale', 'log');
            set(findobj(gca, 'Type', 'Line'), 'LineWidth', 1);
            xlabel('t');
            ylabel('|A|');
            set(fig, 'Position', get(fig, 'Position') + [1 0 0 0] * (get(fig, 'Position') * [0;0;1;0]));
        end
        
        % plots graphiques
        if plotTemporel && verb
            figure;
            for k=1:9
                polarplot([0, meanShape(k)], '-o');
                hold on
            end
            
            title = ['P', num2str(p), 'T', num2str(transient),...
                '_freq=', num2str(meanFreq), '_damp=', num2str(100*zeta)];
            title2 = [title, '_complex'];
            
            fig2 = plotComplexModShape(meanShape, title2);
            
            fig = plotModShape(real(meanShape), title);
            %         fig = plotModShape(imag(meanShape), title);
        end
        
        % progress bar
        if progressBar
            Kprogress = Kprogress + 1;
            waitbar( Kprogress/Nprogress, progressBarFig);
        end
        
        % fin
        while verb
            str = input('continue ? ', 's');
            
            if isequal(str, '0') || isequal(str, 'n') || isequal(str, 'no') || isequal(str, 'non')
                close all
                return
                
            elseif isequal(str, 'test')
                % test
                [t, X] = getData(p, 0);
                boundsT = TransientsTimesP(:, kridge);
                figure;
                plts = plot(t, X(1:9,:));
                plts = transpose(plts);
                
                WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
                    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges, 'CtEdgeEffects', ct,...
                    'XLim', boundsT);
            else
                disp(' ');
                close all
                break
            end
        end
    end
    
end

% progress bar
if progressBar
    close(progressBarFig);
end


if saveFiles
    directory = 'mur silvia\modesOndelette\';
    
    AllModalQuantities = struct('freqs', 0, 'shapes', 0, 'damps', 0);
    
    AllModalQuantities.freqs = allFreqs;
    AllModalQuantities.shapes = allShapes;
    AllModalQuantities.damps = allDamps;
    
    save([directory, 'AllModalQuantities.mat'], 'AllModalQuantities');
end


