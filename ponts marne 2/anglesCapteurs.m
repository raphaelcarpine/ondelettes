clear all

removeNan = false;


%% data

dataFolder = 'C:\Users\carpine\Documents\projets\ponts marne\reprise operations 2021\donnees'; % dossier où les fichier .csv sont
dataFileNames = dir(fullfile(dataFolder, '*niveau*.mat'));
dataFileNames = {dataFileNames.name};

niveauMax = 0;

for kfile = 1:length(dataFileNames)
    
    dataFileName = dataFileNames{kfile};
    
    load(fullfile(dataFolder, dataFileName));
    
%     disp(startDate);
    
    X = X.';
    T = T.';
    
    
    %% tri channels par ordre croissant
    
    [channelNames, I] = sort(channelNames);
    X = X(I, :);
    
    %% remove channels, complete nan
    
    if false
        I = [10, 11];
        %     I = [1:10, 12:14];
        X(I, :) = 0;
    end
    
    if false
        for kc = 1:size(X, 1)
            X(kc, isnan(X(kc, :))) = mean(X(kc, :), 'omitnan');
        end
    end
    
    %% enlever donnees redondantes (double t à cause du mode transmit/log)
    
    
    if 1
        dt0 = mean(diff(T));
        T0 = T;
        X0 = X;
        T = nan(size(T));
        X = nan(size(X));
        
        kt = 1;
        kt0 = 1;
        while kt0 <= length(T0)-1
            if abs(T0(kt0) - T0(kt0+1)) < 1e-4 * dt0
                if all(~isnan(X0(:, kt0)))
                    X(:, kt) = X0(:, kt0);
                elseif all(~isnan(X0(:, kt0+1)))
                    X(:, kt) = X0(:, kt0+1);
                else
                    for kc = 1:size(X0, 1)
                        if ~isnan(X0(kc, kt0))
                            X(kc, kt) = X0(kc, kt0);
                        else
                            X(kc, kt) = X0(kc, kt0+1);
                        end
                    end
                    %                 if any(~isnan(X0(:, kt0:kt0+1)), 'all') % debug
                    %                     disp(X0(:, kt0:kt0+1));
                    %                 end
                    %                 if any(~isnan(X0(:, kt0:kt0+1)), 'all') && any(isnan(X0(:, kt0)) & isnan(X0(:, kt0+1))) % debug
                    %                     disp(X0(:, kt0:kt0+1));
                    %                 end
                end
                T(kt) = T0(kt0);
                kt0 = kt0+2;
            else
                X(:, kt) = X0(:, kt0);
                T(kt) = T0(kt0);
                kt0 = kt0+1;
            end
            kt = kt+1;
        end
        T = T(1:kt-1);
        X = X(:, 1:kt-1);
    end
    
    
    %% enlever nan debut fin
    
    if removeNan
        Nt0 = floor(size(X, 2)/2);
        for Nt1 = Nt0:-1:1
            if any(isnan(X(:, Nt1)))
                Nt1 = Nt1 + 1;
                break
            end
        end
        for Nt2 = Nt0:size(X, 2)
            if any(isnan(X(:, Nt2)))
                Nt2 = Nt2 - 1;
                break
            end
        end
        fprintf('NaN debut : %.2fs\nNaN fin : %.2fs\n', [T(Nt1) - T(1), T(end) - T(Nt2)]);
        X = X(:, Nt1:Nt2);
        T = T(:, Nt1:Nt2);
        
        T = T - T(1);
    end
    
    %% calcul moyenne
    
    X = mean(X, 2, 'omitnan');
    x = X(1:3:end);
    y = X(2:3:end);
    z = X(3:3:end);
    angleSensors = 180/pi * atan(sqrt(x.^2 + y.^2) ./ abs(z));
    
    disp(dataFileName);
    for ks = 1:length(angleSensors)
        disp([channelNames{3*ks-2}(1:end-4), sprintf(': %.2f°', angleSensors(ks))]);
    end
    fprintf('=> max: %.2f°\n\n', max(angleSensors));
    
    niveauMax = max(niveauMax, max(angleSensors));
    
    
    %% deformee
    
    dimensionsShapes2;
    
    shapePlotBridge([angleSensors; nan], dataFileName);
    
end

fprintf('==> global max: %.2f° (%.2f%% distorsion)\n',...
    [niveauMax, 100*(1-cos(pi/180*niveauMax))]);


