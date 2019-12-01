folderName = 'maquetteTMD/donneesTP';

try
    delete([folderName, '/mData.mat']);
catch
end



files = dir(folderName);
filesNames = cell(size(files));
for ind = 1:length(files)
    filesNames{ind} = files(ind).name;
end

firstVar = true;
for ind = 1:length(filesNames)
    fileName = filesNames{ind};
    
    mName = fileName;
    if length(mName) >= 4 && isequal(mName(end-3:end), '.xls')
        mName = mName(1:end-4);
    else
        continue
    end
    
    try
        eval(['data = xlsread(''', fileName, ''');']);
        
        while any(isnan(data(1,:)))
            data = data(2:end,:);
        end
        
        while any(isnan(data(end,:)))
            data = data(1:end-1,:);
        end
        
        
        eval([mName, ' = data;']);        
        
        if firstVar
            save([folderName, '/mData'], mName);
            firstVar = false;
        else
            save([folderName, '/mData'], mName, '-append');
        end
    catch
        warning('problem with variable name');
        continue
    end
end