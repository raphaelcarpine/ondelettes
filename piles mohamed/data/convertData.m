clear all

dataFolder = 'C:\Users\carpine\Documents\projets\experiences affouillement mohamed\data'; % dossier où les fichier .csv sont
saveFolderRaw = fullfile(dataFolder, 'raw mat files'); % dossier où les fichier .mat vont être enregistrés
saveFolderProcessed = dataFolder; % dossier où les fichier .mat vont être enregistrés
extention = '.csv'; % format de l'extension, peut être remplacé par '' si pas d'extension

%% conversion & enregistrement

InfosEssais.pile = [];
InfosEssais.g = [];
InfosEssais.profondeur = [];
InfosEssais.direction = '';
InfosEssais.nbChocs = [];
InfosEssais.resynchronise = false(1, 0);
InfosEssais.timeError = false(1, 0);

for pile = 1:10
    for G = [2 8]
        % recherche des fichiers du dossier
        dataFolder2 = fullfile(dataFolder, sprintf('PILE %u', pile), sprintf('%ug', G));
        files = dir(dataFolder2);
        filesNames = cell(size(files));
        for ind = 1:length(files)
            filesNames{ind} = files(ind).name;
        end
        
        % enregistrement des .csv
        for ind = 1:length(filesNames)
            fileName = filesNames{ind};
            
            mName = fileName;
            if length(mName) >= length(extention) && isequal(mName(end+1-length(extention):end), extention)  % reconnaissance des 'name'.csv
                mName = mName(1:end-length(extention));
            else
                continue
            end
            
            mName = [sprintf('pile%u_', pile), mName(1:end-5), sprintf('_%ug', G)];
            
            % verif fichiers déjà convertis
            files2 = dir(saveFolderProcessed);
            filesNames2 = cell(size(files2));
            for ind = 1:length(files2)
                filesNames2{ind} = files2(ind).name;
            end
            alreadyConverted = false;
            for ind2 = 1:length(filesNames2)
                if isequal(filesNames2{ind2}, [mName, '.mat'])
                    alreadyConverted = true;
                    break
                end
            end
            if alreadyConverted
                continue
            end
            
            % conversion
            
            disp(mName);
            
            D = readtable(fullfile(dataFolder2, fileName), 'VariableNamingRule', 'preserve'); % conversion de 'name'.csv en tableau matlab 'X'
            X = table2array(D(:, 2:end));
            T = table2array(D(:, 1));
            T = T - T(1);
            T = seconds(T);
            channelNames = D.Properties.VariableNames(2:end);
            startDate = D.Time(1);
            startDate.Year = startDate.Year + 2000;
            startDate.TimeZone = 'UTC';
            startDate.TimeZone = 'Europe/Paris';
            
            X = X.';
            T = T.';
            
            % enregistrement raw
            save(fullfile(saveFolderRaw, mName), 'X', 'T', 'channelNames', 'startDate');
            
%             % mise en forme donnees
%             [X, T] = removeRedundantData(X, T);
%             [X, T] = removeNanSignal(X, T);
%             [X, channelNames] = channelNamesConversion(X, channelNames);
%             [~, maxDesync, ~, ~, maxDesyncAfterResync] = testSynchro(T, X, channelNames);
%             fprintf(1 + (maxDesync > 3), 'max desynchronisation: %u*Dt', maxDesync);
%             fprintf(1 + (maxDesyncAfterResync > 3), ' (%u*Dt after re-synchronisation)\n', maxDesyncAfterResync);
%             save(fullfile(saveFolderProcessed, mName), 'X', 'T', 'channelNames', 'startDate', 'maxDesync', 'maxDesyncAfterResync');
            
            % mise en forme donnees
            try
                [X, T] = removeRedundantData(X, T);
                [X, channelNames] = channelNamesConversion(X, channelNames);
                [Tshocks, Xshocks, resynched] = resynchChannels(T, X, channelNames);
                fprintf(1 + isempty(Tshocks), 'nb of shocks: %u ', length(Tshocks));
                if resynched
                    fprintf('(resynchronized)');
                end
                
                timeError = false;
            catch ME
                warning(getReport(ME, 'extended', 'hyperlinks', 'on' ));
                timeError = true;
                Tshocks = {};
                Xshocks = {};
                resynched = false;
            end
            
            fprintf('\n');
            
            save(fullfile(saveFolderProcessed, mName), 'Tshocks', 'Xshocks',...
                'channelNames', 'startDate', 'resynched', 'timeError');
            
            disp(' ');
            
            % infos
            InfosEssais.pile(end+1) = pile;
            InfosEssais.g(end+1) = G;
            InfosEssais.profondeur(end+1) = str2double(mName(7+(pile==10):end-7));
            InfosEssais.direction(end+1) = mName(end-3);
            InfosEssais.nbChocs(end+1) = length(Tshocks);
            InfosEssais.resynchronise(end+1) = resynched;
            InfosEssais.timeError(end+1) = timeError;
        end
    end
end

InfosEssais = structfun(@transpose, InfosEssais, 'UniformOutput', false);
InfosEssais = struct2table(InfosEssais);

writetable(InfosEssais, fullfile(saveFolderProcessed, 'InfosEssais.csv'),'Delimiter',';');

