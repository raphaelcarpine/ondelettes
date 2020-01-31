folderName = 'maquetteTMDesiee/donnees';

try
    delete([folderName, '/mData.mat']); % écrasement de l'ancien fichier .mat
catch
end



files = dir(folderName); % récupération des fichiers du dossier

filesNames = cell(size(files)); % récupération des noms des fichiers
for ind = 1:length(files)
    filesNames{ind} = files(ind).name;
end

firstVar = true;
for ind = 1:length(filesNames)
    fileName = filesNames{ind};
    
    mName = fileName; % nom de la variable
    if length(mName) >= 4 && isequal(mName(end-3:end), '.xls')
        mName = mName(1:end-4);
    else
        continue
    end
    
    try
        eval([mName, '= xlsread(''', fileName, ''');']); % convertion xls/variable matlab
        if firstVar % premier fichier converti
            save([folderName, '/mData'], mName);
            firstVar = false;
        else % suivants
            save([folderName, '/mData'], mName, '-append');
        end
    catch
        warning('problem with variable name');
        continue
    end
end