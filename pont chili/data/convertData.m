% channels :
% accelerometer1 x, y, z (g)
% accelerometer2 x, y, z (g)
% accelerometer3 x, y, z (g)
% displacement * 3*4
% strain + stress * 16



folderName = 'pont chili/data';

% try
%     delete([folderName, '/mData.mat']); % écrasement de l'ancien fichier .mat
% catch
% end

DOFs = 1:53; % all
% DOFs = 1:9; % accelerometers


files = dir(folderName); % récupération des fichiers du dossier

filesNames = cell(size(files)); % récupération des noms des fichiers
for ind = 1:length(files)
    filesNames{ind} = files(ind).name;
end


for ind = 1:length(filesNames)
    fileName = filesNames{ind};
    
    mName = fileName; % nom de la variable
    if length(mName) >= 8 && isequal(mName(end-7:end), '.shm.txt')
        mName = mName(1:end-8);
    else
        continue
    end
    
    try
        file = fopen([folderName, '/', fileName]);
        fileText = fread(file);
        fclose(file);
        
        fileText = char(fileText');
        
        k0 = 1;
        while k0 < length(fileText) && fileText(k0) ~= newline
            k0 = k0+1;
        end
        fileText = fileText(k0+1:end);
        
        fileText = strrep(fileText, '\', '');
        fileText = strrep(fileText, ':', '');
        
        X = str2num(fileText);
        X = transpose(X);
        X = X(DOFs+2, :); % column 1 date, column2 hour
        
        save([folderName, '/', mName], 'X');
    catch
        warning(['problem with conversion; file: ', fileName]);
        continue
    end
end