% affichage et methode de calcul
verb = true;
singleRidgeMode = false;

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


ModesTransients = {{[1; 26], [], [1, 2, 3; 26, 3, 3]},... % P0
    {[1, 2, 3; 3, 1.8, 4.3], [1, 2, 3; 1.8, 1.4, 1.4]},... % P6
    {[1, 2; 17.4, 2.4], [2; 2.5], [1; 17.4]}}; % P7


%%

for ind = 1:3
    p = P(ind);
    freqs = Freqs{ind};
    damps = Damps{ind};
    ModesTransientsP = ModesTransients{ind};
    
    for transient = 1:length(ModesTransientsP)
        modes = ModesTransientsP{transient}(1, :);
        Dfs = ModesTransientsP{transient}(2, :);
        
        for kmode = 1:length(modes)
            mode = modes(kmode);
            
            % donnees mode
            Df = Dfs(kmode);
            f = freqs(mode);
            damp = damps(mode);
            Dt = 1 / (damps * 2*pi*f);
            
            [t, X] = getData(p, transient);
            T = t(end) - t(1);
            
            % choix Q
            [Qmin, Qmax, Qz] = getBoundsQ(f, Df, Dt, T, ct, cf);
            Q = (Qmin + min(Qmax, Qz)) / 2;
            if verb
                disp(['Qmin = ', num2str(Qmin), ' ; Qmax = ', num2str(Qmax), ' ; Qz = ', num2str(Qz)]);
                disp(['Q : ', num2str(Q)]);
                disp(' ');
            end
            
            % calcul modes
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
                
                [time, freq, freqs, shape, amplitude, errors, ridgesNumber] = getModes(ridges, 1);
            end
            
            % calcul moyenne
            
            
        end
    end
end