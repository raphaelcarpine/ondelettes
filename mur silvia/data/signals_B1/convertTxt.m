folderName = 'mur silvia/data/allData'; % chemin du dossier où les fichier .txt sont, et où "mData.mat" va être enregistré
extention = '.txt'; % format de l'extension, peut être remplacé par '' si pas d'extension

%%

% suupression du fichier .mat si il existe déjà (pour reconvertir les
% fichiers si il y a eu un probleme
try
    delete([folderName, '/mData.mat']);
catch
end


% recherche des fichiers du dossier
files = dir(folderName);
filesNames = cell(size(files));
for ind = 1:length(files)
    filesNames{ind} = files(ind).name;
end

% enregistrement des .txt
firstVar = true;
for ind = 1:length(filesNames)
    fileName = filesNames{ind};
    
    mName = fileName;
    if isequal(mName, '.') || isequal(mName, '..') % racines du dossier
        continue
    elseif length(mName) >= length(extention) && isequal(mName(end+1-length(extention):end), extention)  % reconnaissance des 'name'.txt
        mName = mName(1:end-length(extention));
    end
    
    try
        eval([mName, '= dlmread(''', fileName, ''');']); % conversion de 'name'.txt en tableau matlab 'name'
        if firstVar
            save([folderName, '/mData'], mName); % enregistrement
            firstVar = false;
        else
            save([folderName, '/mData'], mName, '-append'); % enregistrement
        end
    catch
        warning('problem with variable name');
        continue
    end
end