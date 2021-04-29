dataFolder = 'C:\Users\carpine\Documents\projets\ponts marne\reprise operations 2021\donnees'; % dossier où les fichier .csv sont
saveFolder = dataFolder; % dossier où les fichier .mat vont être enregistrés
extention = '.csv'; % format de l'extension, peut être remplacé par '' si pas d'extension

%% conversion & enregistrement


% recherche des fichiers du dossier
files = dir(dataFolder);
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
    
    alreadyConverted = false;
    for ind2 = 1:length(filesNames)
        if isequal(filesNames{ind2}, [mName, '.mat'])
            alreadyConverted = true;
            break
        end
    end
    if alreadyConverted
        continue
    end
    
    disp(mName);
    
    try
        D = readtable(fullfile(dataFolder, fileName), 'VariableNamingRule', 'preserve'); % conversion de 'name'.csv en tableau matlab 'X'
        X = table2array(D(:, 2:end));
        T = table2array(D(:, 1));
        T = T - T(1);
        T = seconds(T);
        channelNames = D.Properties.VariableNames(2:end);
        startDate = D.Time(1);
        
        save([saveFolder, '\', mName], 'X', 'T', 'channelNames', 'startDate'); % enregistrement
    catch
        warning('problem with variable name');
        continue
    end
end
