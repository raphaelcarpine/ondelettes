dataFolder = 'piles mohamed\resultats\resultats';
saveFolder = 'piles mohamed\resultats\figures';

saveFigs = 1;
saveGraphsProf = false; % save graphs pour chaque proondeur
pauseBetweenModes = false;

Fs = 4096;

smoothDamping = true;

g = 8;
piles = 10;
directions = 'x';
profondeurs = -15:5:15;

AllAmplFreqPlots = {{}, {}, {}, {}, {}};
AllAmplDampPlots = {{}, {}, {}, {}, {}};
AllProf = {{}, {}, {}, {}, {}};

for pile = piles
    for direction = directions
        for prof = profondeurs
            % file search
            fileNameprefix = sprintf('pile%u_%icm_%c_%ug', [pile, prof, int8(direction), g]);
            files = dir(fullfile(dataFolder, ['*', fileNameprefix, '*.mat']));
            files = {files.name};
            for kf = 1:length(files)
                files{kf} = strsplit(files{kf}, '_');
                files{kf} = files{kf}{5};
                files{kf} = str2double(files{kf}(5:end));
            end
            modes = unique(cell2mat(files));
                
            disp(fileNameprefix);
            
            for mode = modes
                fprintf('mode #%u', mode);
                
                % file search
                fileNameprefix = sprintf('pile%u_%icm_%c_%ug_mode%u', [pile, prof, int8(direction), g, mode]);
                files = dir(fullfile(dataFolder, ['*', fileNameprefix, '*.mat']));
                files = {files.name};
                for kf = 1:length(files)
                    files{kf} = strsplit(files{kf}, '_');
                    files{kf} = files{kf}{6};
                    files{kf} = str2double(files{kf}(6:end-4));
                end
                shocks = unique(cell2mat(files));
                
                AmplFreqPlots = {};
                AmplDampPlots = {};
                ModeShapes = {};
                
                for shock = shocks
                    fileName = sprintf('pile%u_%icm_%c_%ug_mode%u_shock%u.mat',...
                        [pile, prof, int8(direction), g, mode, shock]);
                    load(fullfile(dataFolder, fileName));
                    
                    % freq
                    AmplFreqPlots{end+1} = [abs(amplitudes); freq];
                    
                    % damping
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
                    
                    AmplDampPlots{end+1} = [abs(amplitudes); damp];
                    
                    % mode shapes
                    shape = mean(shapes, 2);
                    if direction == 'x' % orientation similaire
                        shape = shape * (2*(shape(1)>0) - 1);
                    elseif direction == 'y'
                        shape = shape * (2*(shape(2)>0) - 1);
                    end
                    ModeShapes{end+1} = shape;
                end
                
                % 
                AllAmplFreqPlots{mode}{end+1} = AmplFreqPlots;
                AllAmplDampPlots{mode}{end+1} = AmplDampPlots;
                AllProf{mode}{end+1} = prof;
                
                % plots
                figsName = sprintf('pile%u_%icm_%c_%ug_mode%u',...
                        [pile, prof, int8(direction), g, mode]);
                figAD = figure('Name', figsName);
                for ksh = 1:length(AmplDampPlots)
                    plot(AmplDampPlots{ksh}(1, :), 100*AmplDampPlots{ksh}(2, :));
                    hold on
                end
                xlabel('Amplitude [m/s²]');
                ylabel('Damping [%]');
                
                figAF = figure('Name', figsName);
                for ksh = 1:length(AmplFreqPlots)
                    plot(AmplFreqPlots{ksh}(1, :), AmplFreqPlots{ksh}(2, :));
                    hold on
                end
                figAF.Position(1) = figAF.Position(1) + figAF.Position(3);
                xlabel('Amplitude [m/s²]');
                ylabel('Frequency [Hz]');
                
                % mode shape
                shape = mean(cell2mat(ModeShapes), 2);
                I = nonPropIndex(shape);
                FigStr = [figsName, sprintf('; I=%.2f%%', I)];
                
                [fctDefModale, fctDefModaleAnimation] = defModales(pile, channelNames, fileName, prof);
                fctDefModale2 = defModales1D(pile, channelNames, direction, fileName, prof);
                figD0 = fctDefModaleAnimation(shape, FigStr);
                figD0.Position(2) = figD0.Position(2) - figD0.Position(4) - 100;
                figD1 = fctDefModale2(real(shape), FigStr);
                figD1.Position(2) = figD1.Position(2) - figD1.Position(4) - 100;
                figD1.Position(1) = figD1.Position(1) + figD1.Position(3);
                
                % save
                if saveGraphsProf
                    figs = [figAF, figAD, figD1]; % , figD0
                    figs_names = {'ampl_freq', 'ampl_damp', 'mode_shape_1D'}; % , 'mode_shape_3D'
                else
                    figs = figD1;
                    figs_names = {'mode_shape_1D'};
                end
                if saveFigs
                    for kfig = 1:length(figs)
                        saveas(figs(kfig), fullfile(saveFolder, [figsName, '_', figs_names{kfig}]), 'fig');
                        set(figs(kfig),'renderer','Painters');
                        saveas(figs(kfig), fullfile(saveFolder, [figsName, '_', figs_names{kfig}]), 'epsc');
                        saveas(figs(kfig), fullfile(saveFolder, [figsName, '_', figs_names{kfig}]), 'png');
                    end
                    fprintf(' (saved)\n');
                else
                    fprintf('\n');
                end
                
                % wait
                if pauseBetweenModes
                    cont = input('continue? ', 's');
                    switch cont
                        case {'', 'y', 'yes', 'oui'}
                            close all
                        otherwise
                            return
                    end
                else
                    close all
                end
            end
            
            fprintf('\n');
        end
    end
end

%% plot all curves

for km = 1:length(AllAmplFreqPlots) % mode
    if isempty(AllAmplFreqPlots{km})
        break
    end
    
    % freq
    fig = figure('Name', sprintf('freq: pile%u_%c_%ug_mode%u', [pile, int8(direction), g, km]));
    for kp = 1:length(AllAmplFreqPlots{km}) % prof
        AmplFreqPlots = AllAmplFreqPlots{km}{kp};
        prof = AllProf{km}{kp};
        
        for ks = 1:length(AmplFreqPlots) % shock
            AmplFreqPlots{ks} = [AmplFreqPlots{ks}, [nan; nan]];
        end
        AmplFreqPlots = cell2mat(AmplFreqPlots);
        plot(AmplFreqPlots(1, :), AmplFreqPlots(2, :), 'DisplayName', sprintf('%i cm', prof));
        xlabel('Amplitude [m/s²]');
        ylabel('Frequency [Hz]');
        hold on
    end
%     plot([10 10], [nan nan], 'HandleVisibility', 'off');
    ax = gca;
    ax.XLim(1) = 0;
    legend();
    
    figName = sprintf('pile%u_%c_%ug_mode%u', [pile, int8(direction), g, km]);
    fig_names = {'ampl_freq', 'ampl_damp'};
    if saveFigs
        saveas(fig, fullfile(saveFolder, [figName, '_', fig_names{1}]), 'fig');
        set(fig,'renderer','Painters');
        saveas(fig, fullfile(saveFolder, [figName, '_', fig_names{1}]), 'epsc');
        saveas(fig, fullfile(saveFolder, [figName, '_', fig_names{1}]), 'png');
        fprintf(' (saved)\n');
    end
    
    % damping
    fig = figure('Name', sprintf('damping: pile%u_%c_%ug_mode%u', [pile, int8(direction), g, km]));
    for kp = 1:length(AllAmplFreqPlots{km}) % prof
        AmplDampPlots = AllAmplDampPlots{km}{kp};
        prof = AllProf{km}{kp};
        
        for ks = 1:length(AmplDampPlots) % shock
            AmplDampPlots{ks} = [AmplDampPlots{ks}, [nan; nan]];
        end
        AmplDampPlots = cell2mat(AmplDampPlots);
        plot(AmplDampPlots(1, :), 100*AmplDampPlots(2, :), 'DisplayName', sprintf('%i cm', prof));
        xlabel('Amplitude [m/s²]');
        ylabel('Damping [%]');
        hold on
    end
%     plot([10 10], [100 110], 'HandleVisibility', 'off', 'Visible', 'off');
    ax = gca;
    ax.XLim(1) = 0;
    ax.YLim(1) = 0;
    legend();
    
    if saveFigs
        saveas(fig, fullfile(saveFolder, [figName, '_', fig_names{2}]), 'fig');
        set(fig,'renderer','Painters');
        saveas(fig, fullfile(saveFolder, [figName, '_', fig_names{2}]), 'epsc');
        saveas(fig, fullfile(saveFolder, [figName, '_', fig_names{2}]), 'png');
        fprintf(' (saved)\n');
    end
    
end