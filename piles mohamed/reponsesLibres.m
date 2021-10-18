dataFolder  = 'C:\Users\carpine\Documents\projets\experiences affouillement mohamed\data';
saveFolder = 'piles mohamed\resultats\resultats';
Fs = 4096;

saveResult = 0;
plotResult = 1;
manualMode = 1;
saveResult = saveResult && ~manualMode;
plotResult = plotResult && ~manualMode;

selectChannelsDir = true;
MotherWavelet = 'cauchy'; % 'morlet', 'cauchy'
MultiSignalModeAbsValue = true;
MultiSignalModeAbsValueModeAmpl = true;
smoothDamping = true;
MaxSlopeRidge = 2;

g = 8;

piles = 3;
directions = 'x';
profondeurs = 5:5:30;

for pile = piles
    for direction = directions
        for prof = profondeurs
            fileName = sprintf('pile%u_%icm_%c_%ug', [pile, prof, int8(direction), g]);
            load(fullfile(dataFolder, fileName));
            fprintf(fileName);
            fprintf(' (%u shocks)\n', length(Tshocks));
            
            if selectChannelsDir
                for ks = 1:length(Xshocks)
                    if ks < length(Xshocks)
                        Xshocks{ks} = selectChannelsDirection(Xshocks{ks}, channelNames, direction);
                    elseif ks == length(Xshocks)
                        [Xshocks{ks}, channelNames] = selectChannelsDirection(Xshocks{ks}, channelNames, direction);
                    end
                end
            end
            
            for shock = 1:length(Tshocks)
                fprintf('shock #%u\n', shock);
                
                [freqsModes, dampingsModes, other_freqs0, lag] = getModalParam(pile, prof, direction, g);
                modes = 1:length(freqsModes);
                freqsModes2 = [0, freqsModes, inf];
                
                T = Tshocks{shock};
                X = Xshocks{shock};
                T = T - T(1);
                if ~isempty(lag)
                    X = X(:, T >= lag);
                    T = T(T >= lag);
                    T = T - T(1);
                end
                X = X - mean(X, 2);
                
                for mode = modes
                    otherFreqs = [cell2mat(freqsModes([1:mode-1, mode+1:end])), other_freqs0, inf]; %-freqsModes{1}(1)
                    otherFreqs1 = [0, otherFreqs(otherFreqs < freqsModes{mode}(1))];
                    otherFreqs2 = otherFreqs(otherFreqs > freqsModes{mode}(end));
                    fmin = max((max(otherFreqs1) + freqsModes{mode}(1))/2, freqsModes{mode}(1) - 5);
                    fmax = min((min(otherFreqs2) + freqsModes{mode}(end))/2, freqsModes{mode}(end) + 5);
                    Df = min([abs(otherFreqs - freqsModes{mode}(1)),...
                        abs(otherFreqs - freqsModes{mode}(end))]);
                    [Qmin, Qmax, Qz] = getBoundsQ(mean(freqsModes{mode}), Df,...
                        1/(mean(dampingsModes{mode})*2*pi*mean(freqsModes{mode})), inf, 3, 5);
                    if Qmin > Qz
                        warning(sprintf('%.1f = Qmin > Qz = %.1f', [Qmin, Qz]));
                    end
                    Q = max(Qmin, 2);
                    NbFreq = 300;
                    
                    fprintf('    mode #%u (Q=%.2f)', [mode, Q]);
                    
                    % mode shapes plots
                    [fctDefModale, fctDefModaleAnimation] = defModales(pile, channelNames, fileName, prof);
                    fctDefModale2 = defModales1D(pile, channelNames, direction, fileName, prof);
                    
                    % wavelet plot
                    if manualMode
                        fig = figure('Name', fileName);
                        plt = plot(T, X);
                        legend(channelNames);
                        selectLine();
                        
                        WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'RemoveMean', true,...
                            'RealShapePlot', fctDefModale2, 'FourierScale', 'lin', 'MotherWavelet', MotherWavelet,...
                            'AnimatedShapePlot', fctDefModaleAnimation, 'StopRidgeWhenIncreasing', true,...
                            'MultiSignalMode', true, 'MultiSignalModeAbsValue', MultiSignalModeAbsValue,...
                            'MaxSlopeRidge', MaxSlopeRidge);
                        waitfor(fig);
                    end
                    
                    % ridge extraction
                    [t, freq, shapes, amplitudes] = ...
                        getModesSingleRidge(T, X, Q, fmin, fmax, NbFreq,...
                        'NbMaxParallelRidges', inf, 'NbMaxRidges', 1, 'MinModu', 0,...
                        'StopWhenIncreasing', true, 'MaxSlopeRidge', MaxSlopeRidge,...
                        'MotherWavelet', MotherWavelet, 'MultiSignalModeAbsValue', MultiSignalModeAbsValue,...
                        'MultiSignalModeAbsValueModeAmpl', MultiSignalModeAbsValueModeAmpl);
                    t = t{1};
                    freq = freq{1};
                    shapes = shapes{1};
                    amplitudes = amplitudes{1};
                    
                    % plot
                    if plotResult
                        FigStr = sprintf('pile%u; %icm; %c; %ug; mode%u; shock%u',...
                            [pile, prof, int8(direction), g, mode, shock]);
                        
                        % CWT
                        WvltFreq = linspace(fmin, fmax, NbFreq);
                        CWT = zeros(NbFreq, size(X, 2));
                        for kch = 1:size(X, 1)
                            if MultiSignalModeAbsValue
                                CWT = CWT +  abs(WvltComp(T, X(kch, :), WvltFreq, Q, 'MotherWavelet', MotherWavelet)).^2;
                            else
                                CWT = CWT +  WvltComp(T, X(kch, :), WvltFreq, Q, 'MotherWavelet', MotherWavelet).^2;
                            end
                        end
                        WvltPlot2(T, WvltFreq, CWT, 'abs', Q, 3, MotherWavelet, 'log', 'sum CWT', '', 'lin');
                        
                        % ridge
                        figure('Name', FigStr);
                        plot(t, freq);
                        xlabel('Time [s]');
                        ylabel('Frequency [Hz]');
                        
                        figure('Name', FigStr);
                        plot(t, abs(amplitudes));
                        set(gca, 'YScale', 'log');
                        xlabel('Time [s]');
                        ylabel('Amplitude [m/s²]');
                        
                        figure('Name', FigStr);
                        plot(abs(amplitudes), freq);
                        xlabel('Amplitude [m/s²]');
                        ylabel('Frequency [Hz]');
                        
                        % plot damping
                        damp0 = -diff(log(abs(amplitudes)))*Fs;
                        damp0 = ([damp0(1), damp0] + [damp0, damp0(end)])/2;
                        if smoothDamping
                            deltaT = Q / (2*pi*mean(freq));
                            deltaT = deltaT/5;
                            Nhalf = ceil(5*deltaT*Fs);
                            Nhalf = min(Nhalf, length(damp0));
                            filterCoeffs = nan(1, 2*Nhalf+1);
                            for kfilter = 1:length(filterCoeffs)
                                filterCoeffs(kfilter) = exp( -((kfilter-1-Nhalf)/Fs)^2 / (2*deltaT^2));
                            end
                            filterCoeffs = filterCoeffs / sum(filterCoeffs);
                            damp0 = [damp0(1)*ones(1, Nhalf), damp0, damp0(end)*ones(1, Nhalf)];
                            damp0 = filter(filterCoeffs, 1, damp0);
                            damp0 = damp0(2*Nhalf+1:end);
                        end
                        damp0 = damp0/(2*pi);
                        damp = damp0 ./ sqrt(damp0.^2 + freq.^2);
                        fig = figure('Name', FigStr);
                        fig.Position(1) = fig.Position(1) + fig.Position(3);
                        plot(abs(amplitudes), 100*damp);
                        xlabel('Amplitude [m/s²]');
                        ylabel('Damping [%]');
                        
                        % mode shapes
                        shape = mean(shapes, 2);
                        I = nonPropIndex(shape);
                        FigStr = [FigStr, sprintf('; I=%.2f%%', I)];
                        fig = fctDefModale2(real(shape), FigStr);
                        fig.Position(2) = fig.Position(2) - fig.Position(4) - 80;
                        fig = fctDefModaleAnimation(shape, FigStr);
                        fig.Position(2) = fig.Position(2) - fig.Position(4) - 80;
                        fig.Position(1) = fig.Position(1) + fig.Position(3);
                    end
                    
                    % saving
                    if saveResult
                        fileName = sprintf('pile%u_%icm_%c_%ug_mode%u_shock%u',...
                            [pile, prof, int8(direction), g, mode, shock]);
                        save(fullfile(saveFolder, fileName), 't', 'freq', 'shapes', 'amplitudes',...
                            'Q', 'channelNames', 'selectChannelsDir', 'MotherWavelet', 'MultiSignalModeAbsValue');
                        fprintf(' (saved)\n');
                    else
                        fprintf('\n');
                    end
                    
                    % wait
                    if plotResult
                        cont = input('continue? ', 's');
                        switch cont
                            case {'', 'y', 'yes', 'oui'}
                                close all
                            otherwise
                                return
                        end
                    end
                end
            end
            
            fprintf('\n');
        end
    end
end