load('mur silvia\modesFourier\ModesLMS.mat');

ModesWvlt = struct([]);

%%
% affichage et methode de calcul
verb = true;
singleRidgeMode = false;
plotTemporel = true;
% test = false;
save = false;

% choix precision
ct = 3;
cf = 5;

%% data

P = [0, 6, 7];

Freqs = {[8.35, 33.95, 36.77],...
    [11.12, 32.79, 37.76],...
    [10.97, 28.34, 34.22]};

Damps = {[2.16, 0.47, 0.44] * 0.01,...
    [0.57, 0.72, 0.99] * 0.01,...
    [1.03, 0.68, 0.71] * 0.01};


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
TransientsTimes{1} = [TransientsTimes{1}, [1267.35, 1267.35; 1269.95, 1269.95]];
TransientsDeltaF{1} = [TransientsDeltaF{1}, 2.9, 2.9];

%P0T2
TransientsModes{1} = [TransientsModes{1}, 1, 2, 3];
TransientsNumbers{1} = [TransientsNumbers{1}, 2, 2, 2];
TransientsTimes{1} = [TransientsTimes{1}, [1450.2, 1450.2, 1450.2; 1452.5, 1453, 1452.2]];
TransientsDeltaF{1} = [TransientsDeltaF{1}, 2.7, 2.9, 2.9];

%P6T1
TransientsModes{2} = [TransientsModes{2}, 1, 2, 3];
TransientsNumbers{2} = [TransientsNumbers{2}, 1, 1, 1];
TransientsTimes{2} = [TransientsTimes{2}, [238, 238, 238; 243, 243, 239.5]];
TransientsDeltaF{2} = [TransientsDeltaF{2}, 4.7, 1.5, 4.9];

%P6T2
TransientsModes{2} = [TransientsModes{2}, 2];
TransientsNumbers{2} = [TransientsNumbers{2}, 2];
TransientsTimes{2} = [TransientsTimes{2}, [444; 447]];
TransientsDeltaF{2} = [TransientsDeltaF{2}, 1.4];



%%

for indp = 1:3
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
        Q = (Qmin + min(Qmax, Qz)) / 2;
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
        
        % comparaison fourier
        errorFreq = abs(meanFreq - ModesLMS(indp, mode).freq) / ModesLMS(indp, mode).freq;
        
        shapeF = ModesLMS(indp, mode).shape;
        errorShape = meanShape - shapeF;
        errorShape = sqrt (sum((abs(errorShape).^2)) / sum((abs(shapeF).^2)));
        
        if verb
            disp(['freq : ', num2str(meanFreq), ' ; error : ', num2str(100*errorFreq), '%']);
            disp(['shape error : ', num2str(100*errorShape), '%']);
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
        figure;
        for k=1:9
            polarplot([0, meanShape(k)], '-o');
            hold on
        end
        
        title = ['P', num2str(p), 'T', num2str(transient),...
            '_freq=', num2str(meanFreq), '_damp=', num2str(zeta)];
        title2 = [title, '_complex'];
        
        fig2 = plotComplexModShape(meanShape, title2);
        
        fig = plotModShape(real(meanShape), title);
        
        % enregistrement
        if save
            directory = 'mur silvia\modesOndelette\';
            savefig(fig, [directory, 'save\', title, '.fig']);
            saveas(fig, [directory, 'save\', title, '.png']);
            savefig(fig2, [directory, 'save\', title2, '.fig']);
            saveas(fig2, [directory, 'save\', title2, '.png']);
        end
        
        % fin
        while true
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
                close all
                break
            end
        end
        
    end
end