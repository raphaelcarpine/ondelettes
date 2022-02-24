folder = 'C:\Users\carpine\Documents\projets\simulations elements finis non lin\data';

files = dir([folder, '/**/*simul*.mat']);
fileNames = {files.name};
fileFolders = {files.folder};

nbs0 = flip(12:26);
nbs1 = flip([(12:21)+5, (22:26)+10]);

for kn = 1:length(nbs0)
    nb0 = nbs0(kn);
    nb1 = nbs1(kn);
    for kf = 1:length(fileNames)
        fileFolder = fileFolders{kf};
        fileName = fileNames{kf};
        if strcmp(fileName(1:7), sprintf('simul%d', nb0))
            fileName2 = fileName;
            fileName2(1:7) = sprintf('simul%d', nb1);
            movefile(fullfile(fileFolder, fileName), fullfile(fileFolder, fileName2));
        end
    end
end