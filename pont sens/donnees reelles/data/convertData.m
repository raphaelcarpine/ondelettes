dataFolder = 'C:\Users\carpine\Documents\projets\article bruno\Pont de Sens (données)\donnees_brutes\24_06_03_Avant resserrage\données'; % dossier où les fichier .txt sont
saveFolder = 'pont sens\donnees reelles\data'; % dossier où les fichier .mat vont être enregistrés
extention = '.TXT'; % format de l'extension, peut être remplacé par '' si pas d'extension

%% conversion & enregistrement


% recherche des fichiers du dossier
files = dir(dataFolder);
filesNames = cell(size(files));
for ind = 1:length(files)
    filesNames{ind} = files(ind).name;
end

% enregistrement des .txt
for ind = 1:length(filesNames)
    fileName = filesNames{ind};
    
    mName = fileName;
    if length(mName) >= length(extention) && isequal(mName(end+1-length(extention):end), extention)  % reconnaissance des 'name'.txt
        mName = mName(1:end-length(extention));
    else
        continue
    end
    
    try
        eval([mName, '= dlmread(''', dataFolder, '\' fileName, ''');']); % conversion de 'name'.txt en tableau matlab 'name'
        save([saveFolder, '\', mName], mName); % enregistrement
    catch
        warning('problem with variable name');
        continue
    end
end