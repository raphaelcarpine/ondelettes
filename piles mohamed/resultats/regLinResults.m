dataFolder = 'piles mohamed\resultats\resultats';
saveFolder = 'piles mohamed\resultats\reg lin';

plotRegs = false;
saveRegs = true;

Fs = 4096;

smoothDamping = true;

g = 8;
pile = 10;
direction = 'x';

Profondeurs = -15:5:15;
Profondeurs = Profondeurs.';

% table
Frequence_0 = nan(size(Profondeurs));
Frequence_coeff = nan(size(Profondeurs));
Amortissement_0 = nan(size(Profondeurs));
Amortissement_coeff = nan(size(Profondeurs));
regLinTable = table(Profondeurs, Frequence_0, Frequence_coeff, Amortissement_0, Amortissement_coeff);
regLinTables = {regLinTable, regLinTable, regLinTable, regLinTable}; % mode 1n mode 2, mode 3, mode 4,...

Nmodes = 0;

for kprof = 1:length(Profondeurs)
    prof = Profondeurs(kprof);
    
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
    modes = reshape(modes, [1, length(modes)]);
    
    disp(fileNameprefix);
    
    for mode = modes
        fprintf('mode #%u\n', mode);
        Nmodes = max(Nmodes, mode);
        
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
        shocks = reshape(shocks, [1, length(shocks)]);
        
        AmplFreqPlots = {};
        AmplDampPlots = {};
        
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
        end
        
        % linear reg
        [regLinTables{mode}.Frequence_0(kprof), regLinTables{mode}.Frequence_coeff(kprof)] = regLinFunc(AmplFreqPlots);
        [regLinTables{mode}.Amortissement_0(kprof), regLinTables{mode}.Amortissement_coeff(kprof)] = regLinFunc(AmplDampPlots);
        
        % plot
        if plotRegs
            figsName = sprintf('pile%u_%icm_%c_%ug_mode%u',...
                [pile, prof, int8(direction), g, mode]);
            figAD = figure('Name', figsName);
            for ksh = 1:length(AmplDampPlots)
                plot(AmplDampPlots{ksh}(1, :), 100*AmplDampPlots{ksh}(2, :));
                hold on
            end
            x0 = get(gca, 'XLim');
            plot(x0, 100*(regLinTables{mode}.Amortissement_0(kprof) + regLinTables{mode}.Amortissement_coeff(kprof) * x0), 'r--');
            xlabel('Amplitude [m/s²]');
            ylabel('Damping [%]');
            
            figAF = figure('Name', figsName);
            for ksh = 1:length(AmplFreqPlots)
                plot(AmplFreqPlots{ksh}(1, :), AmplFreqPlots{ksh}(2, :));
                hold on
            end
            x0 = get(gca, 'XLim');
            plot(x0, regLinTables{mode}.Frequence_0(kprof) + regLinTables{mode}.Frequence_coeff(kprof) * x0, 'r--');
            xlabel('Amplitude [m/s²]');
            ylabel('Frequency [Hz]');
            
            % wait
            cont = input('continue? ', 's');
            switch cont
                case {'', 'y', 'yes', 'oui'}
                    close all
                otherwise
                    return
            end
        end
    end
    
    fprintf('\n');
end

regLinTables = regLinTables(1:Nmodes);

for mode = 1:Nmodes
    fprintf('\nmode #%u\n', mode);
    regLinTable2 = regLinTables{mode};
    regLinTable2.Amortissement_0 = 100*regLinTable2.Amortissement_0;
    regLinTable2.Amortissement_coeff = 100*regLinTable2.Amortissement_coeff;
    
    disp(regLinTable2);
end

if saveRegs
    for mode = 1:Nmodes
        fileName = sprintf('pile%u_%c_%ug_mode%u_reglin.xlsx', [pile, int8(direction), g, mode]);
        writetable(regLinTables{mode}, fullfile(saveFolder, fileName));
        disp('saved');
    end
end