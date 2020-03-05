load('mur silvia\modesFourier\ModesLMS.mat');

ModesWvlt = struct([]);

%%
% affichage et methode de calcul
verb1 = true; % résultats globaux
verb2 = false; % résultats transients
singleRidgeMode = true;
plotTemporel = false;
plotGlobal = true;
save = false;

% choix precision
ct = 3;
cf = 5;

%% moyennes et ecarts types

% nombre de modes pour P0, P6 et P7
nbModes = [3, 4, 3];

% moyennes et ecarts types
for indp = 1:3
    for mode = 1:nbModes(indp)
        allFreqs{indp}{mode} = [];
        allShapes{indp}{mode} = [];
        allDamps{indp}{mode} = [];
    end
end

meanFreqs = {[], [], []};
meanShapes = {[], [], []};
meanDamps = {[], [], []};
nbTransients = {[], [], []};
stdFreqs = {[], [], []};
stdShapes = {[], [], []};
stdDamps = {[], [], []};

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


%P0T1
TransientsModes{1} = [TransientsModes{1}, 2, 3];
TransientsNumbers{1} = [TransientsNumbers{1}, 1, 1];
TransientsTimes{1} = [TransientsTimes{1}, [1218.1, 1218.1; 1223, 1219.4]];
TransientsDeltaF{1} = [TransientsDeltaF{1}, 2.8, 2.8];

%P0T2
TransientsModes{1} = [TransientsModes{1}, 2, 3];
TransientsNumbers{1} = [TransientsNumbers{1}, 2, 2];
TransientsTimes{1} = [TransientsTimes{1}, [1267.35, 1267.35; 1269.95, 1269.95]];
TransientsDeltaF{1} = [TransientsDeltaF{1}, 2.8, 2.8];

%P0T3
TransientsModes{1} = [TransientsModes{1}, 1, 2, 3];
TransientsNumbers{1} = [TransientsNumbers{1}, 3, 3, 3];
TransientsTimes{1} = [TransientsTimes{1}, [1450.2, 1450.2, 1450.2; 1452.5, 1453, 1452.2]];
TransientsDeltaF{1} = [TransientsDeltaF{1}, 2.7, 2.8, 2.8];

%P6T1
TransientsModes{2} = [TransientsModes{2}, 1, 2, 4];
TransientsNumbers{2} = [TransientsNumbers{2}, 1, 1, 1];
TransientsTimes{2} = [TransientsTimes{2}, [238, 238, 238; 243, 241, 239.5]];
TransientsDeltaF{2} = [TransientsDeltaF{2}, 4.7, 1.7, 4.9];

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
TransientsTimes{3} = [TransientsTimes{3}, [286.7, 286.7; 288, 288]];
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

if save
    delete('mur silvia\modesOndelette\save\*');
end

for indp = 1:3
    p = P(indp);
    freqs = Freqs{indp};
    damps = Damps{indp};
    
    TransientsModesP = TransientsModes{indp};
    TransientsNumbersP = TransientsNumbers{indp};
    TransientsTimesP = TransientsTimes{indp};
    TransientsDeltaFP = TransientsDeltaF{indp};
    
    if verb2
        disp(' ');
        disp(['~~~~~~ P', num2str(p), ' ~~~~~~']);
        disp(' ');
    end
    
    for kridge = 1:length(TransientsNumbersP)
        transient = TransientsNumbersP(kridge);
        mode = TransientsModesP(kridge);
        
        if verb2
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
        if verb2
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
                    'NbMaxParallelRidges', MaxRidges, 'NbMaxRidges', MaxParallelRidges);
            end
            
            [time, freq, ~, shape, amplitude, errors, ridgesNumber] = getModes(ridges, 1);
        end
        
        % calcul moyenne
        if isempty(time)
            warning(['no ridge, P', num2str(p), 'T', num2str(transient), 'm', num2str(mode)]);
            continue
        end
        
        if length(time) > 1
            warning('multiple ridges');
        end
        
        meanFreq = mean(freq{1});
        shapeNoNan = shape{1};
        shapeNoNan(isnan(shapeNoNan)) = 0;
        meanShape = mean(shapeNoNan, 2);
        
        funReg = @(param) param(1).*exp(-param(2)*(time{1}-time{1}(1))) - abs(amplitude{1});
        options = optimoptions(@lsqnonlin, 'Display', 'off');
        A0lambda = lsqnonlin(funReg, [abs(amplitude{1}(1)), damp], [], [],...
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
        
        if verb2
            disp(['freq : ', num2str(meanFreq), ' ; freq lms : ', num2str(ModesLMS(indp, mode).freq), ' ; error : ', num2str(100*errorFreq), '%']);
            disp(['MAC : ', num2str(100*mac), '%, ', 'shape error : ', num2str(100*errorShape), '%']);
            disp(['amort : ', num2str(100*zeta), '% ; amort lms : ', num2str(100*ModesLMS(indp, mode).damping), '%']);
            disp(['I : ', num2str(100*nonPropIndex(meanShape)), '%',...
                ' ; I lms : ', num2str(100*nonPropIndex(ModesLMS(indp, mode).shape)), '%']);
        end
        
        % plots temporels
        figs = [];
        if plotTemporel
            figure;
            plot(time{1}, (angle(shape{1}*exp(-1i*pi/2)) + pi/2) * 180/pi);
            hold on
            plot(time{1}, zeros(size(time{1})), 'black--');
            plot(time{1}, 180*ones(size(time{1})), 'black--');
            ylim([-90, 270]);
            xlabel('t');
            ylabel('arg(T)');
            
            figure;
            plot(time{1}, real(shape{1}));
            xlabel('t');
            ylabel('Re(T)');
            figure;
            plot(time{1}, imag(shape{1}));
            xlabel('t');
            ylabel('Im(T)');
            
            figure;
            plot(time{1}, freq{1});
            xlabel('t');
            ylabel('f');
            
            figure;
            hold on
            plot(time{1}, abs(amplitude{1}));
            plot(time{1}, A0*exp(-lambda*(time{1}-time{1}(1))), 'r--');
            set(gca, 'YScale', 'log');
            xlabel('t');
            ylabel('|A|');
        end
        
        % plots graphiques
        if plotTemporel
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
        
        % fin
        while verb2
            str = input('continue ? ', 's');
            
            if isequal(str, '0') || isequal(str, 'n') || isequal(str, 'no') || isequal(str, 'non')
                close all
                return
                
            elseif isequal(str, 'test')
                % test
                figure;
                plts = plot(t, X(1:9,:));
                plts = transpose(plts);
                
                WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
                    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges, 'CtEdgeEffects', ct);
                
            else
                disp(' ');
                close all
                break
            end
        end
    end
    
    % moyennes
    if verb1
        disp(' ');
        disp(' ');
        disp(['~~~~~ P', num2str(p), ' : mean, std']);
    end
    
    for mode = 1:nbModes(indp)
        meanFreqs{indp}(mode) = mean(allFreqs{indp}{mode});
        stdFreqs{indp}(mode) = std(allFreqs{indp}{mode});
        meanDamps{indp}(mode) = mean(allDamps{indp}{mode});
        stdDamps{indp}(mode) = std(allDamps{indp}{mode});
        
        meanShapes{indp}(mode, :) = mean(allShapes{indp}{mode}, 1);
        stdShapes{indp}(mode, :) = std( real( allShapes{indp}{mode}), 0, 1) + 1i*std( imag( allShapes{indp}{mode}), 0, 1);
        
        nbTransients{indp}(mode) = length(allFreqs{indp}{mode});
        
        meanFreq = meanFreqs{indp}(mode);
        meanShape = transpose(meanShapes{indp}(mode, :));
        zeta = meanDamps{indp}(mode);
        
        % erreur stat
        errorFreq = stdFreqs{indp}(mode) / sqrt(nbTransients{indp}(mode));
        errorDamp = stdDamps{indp}(mode) / sqrt(nbTransients{indp}(mode));
        errorShape = norm( stdShapes{indp}(mode, :)) / sqrt(nbTransients{indp}(mode));
        errorShape = errorShape / norm(meanShape);
        shapeI = norm( imag( meanShape));
        errorShapeI = norm( imag( stdShapes{indp}(mode, :))) / sqrt(nbTransients{indp}(mode));
        
        % comparaison fourier
        errorFreqF = abs(meanFreq - ModesLMS(indp, mode).freq) / ModesLMS(indp, mode).freq;
        
        shapeF = ModesLMS(indp, mode).shape;
        errorShapeF = meanShape - shapeF;
        errorShapeF = sqrt( errorShapeF'*errorShapeF / (shapeF'*shapeF));
        
        % modal assurance criterion
        mac = abs(meanShape'*shapeF)^2 / ((meanShape'*meanShape) * (shapeF'*shapeF));
        
        if verb1
            disp(' ');
            disp(['~ mode', num2str(mode), ' (', num2str(nbTransients{indp}(mode)), ' transients)']);
            disp(['freq : ', num2str(meanFreq), ' ; error : ', num2str(100*errorFreq/meanFreq), '%',...
                ' ; freq lms : ', num2str(ModesLMS(indp, mode).freq), ' ; error lms : ', num2str(100*errorFreqF), '%']);
            disp(['MAC : ', num2str(100*mac), '% ; ', ' ; shape error : ', num2str(100*errorShape), '%',...
                ' ; shape error lms : ', num2str(100*errorShapeF), '%']);
            disp(['imaginary shape : ', num2str(shapeI), ' ; ', 'imaginary shape error : ',...
                num2str(100*errorShapeI/shapeI), '%']);
            disp(['amort : ', num2str(100*zeta), '% ; error : ', num2str(errorDamp),...
                ' ; amort lms : ', num2str(100*ModesLMS(indp, mode).damping), '%']);
            
            disp(['I : ', num2str(100*nonPropIndex(meanShape)), '%',...
                ' ; I lms : ', num2str(100*nonPropIndex(ModesLMS(indp, mode).shape)), '%']);
        end
        
        if plotGlobal
            
            
            
            % enregistrement
            if save
                directory = 'mur silvia\modesOndelette\save\';
                savefig(fig, [directory, title, '.fig']);
                saveas(fig, [directory, title, '.png']);
                savefig(fig2, [directory, title2, '.fig']);
                saveas(fig2, [directory, title2, '.png']);
            end
        end
        
    end
    
end


