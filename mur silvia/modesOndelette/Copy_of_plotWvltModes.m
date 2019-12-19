% affichage et methode de calcul
verb = true;
singleRidgeMode = false;
plotTemporel = true;
test = false;
save = false;

% choix precision
ct = 3;
cf = 5;

%% data

P = [0, 6, 7];

Freqs = {[8.35, 33.95, 36.77],...
    [11.12, 32.79, 37.76],...
    [10.96, 31.63, 36.98]};

Damps = {[2.16, 0.47, 0.44] * 0.01,...
    [0.57, 0.72, 0.99] * 0.01,...
    [0.52, 1.16, 1.06] * 0.01};


% ModesTransients = {{[1; 26], [], [1, 2, 3; 26, 3, 3]},... % P0
%     {[1, 2, 3; 3, 1.8, 4.3], [1, 2, 3; 1.8, 1.4, 1.4]},... % P6
%     {[1, 2; 17.4, 2.4], [2; 2.5], [1; 17.4]}}; % P7

ModesTransients = {{[1; 26], [], [1, 2, 3; 26, 3, 3]},... % P0
    {[1, 2, 3; 3, 1.8, 4.3], [2; 1.4]},... % P6
    {[2; 2.4], [2; 2.5], [1; 17.4]}}; % P7


%%

for ind = 1:3
    p = P(ind);
    freqs = Freqs{ind};
    damps = Damps{ind};
    ModesTransientsP = ModesTransients{ind};
    
    if verb
        disp(' ');
        disp(['~~~~~~ P', num2str(p), ' ~~~~~~']);
    end
    
    for transient = 1:length(ModesTransientsP)
        
        if verb
            disp(['~~~ T', num2str(transient)]);
        end
        
        for kmode = 1:size(ModesTransientsP{transient}, 2)
            modes = ModesTransientsP{transient}(1, :);
            Dfs = ModesTransientsP{transient}(2, :);
            
            mode = modes(kmode);
            
            % donnees mode
            Df = Dfs(kmode);
            f = freqs(mode);
            damp = damps(mode);
            Dt = 1 / (damp * 2*pi*f);
            
            [t, X] = getData(p, transient);
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
            meanShape = mean(shape{1}, 2);
            
            funReg = @(param) param(1).*exp(-param(2)*(time{1}-time{1}(1))) - abs(amplitude{1});
            options = optimoptions(@lsqnonlin, 'Display', 'off');
            A0lambda = lsqnonlin(funReg, [abs(amplitude{1}(1)), damp], [], [],...
                options);
            A0 = A0lambda(1);
            lambda = A0lambda(2);
            zeta = lambda / (2*pi*meanFreq);
            
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
end