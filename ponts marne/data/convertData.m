folderName = 'ponts marne/data/tests';

% try
%     delete([folderName, '/mData.mat']); % écrasement de l'ancien fichier .mat
% catch
% end

DOFs = 1:6;


files = dir(folderName); % récupération des fichiers du dossier

filesNames = cell(size(files)); % récupération des noms des fichiers
for ind = 1:length(files)
    filesNames{ind} = files(ind).name;
end


for ind = 1:length(filesNames)
    fileName = filesNames{ind};
    
    mName = fileName; % nom de la variable
    if length(mName) >= 4 && isequal(mName(end-3:end), '.ASC')
        mName = mName(1:end-4);
    else
        continue
    end
    
    try
        file = fopen([folderName, '/', fileName]);
        fileText = fread(file);
        fclose(file);
        
        fileText = char(fileText');
        
        fileText = strrep(fileText, ',', '.');
        fileText = strrep(fileText, ':', '');
        fileText = strrep(fileText, '/', '');
        fileText = strrep(fileText, ' ', '');
        
        X = str2num(fileText);
        X = transpose(X);
        X = X(DOFs+1, :);
        
        save([folderName, '/', mName], 'X');
    catch
        warning(['problem with conversion; file: ', fileName]);
        continue
    end
end