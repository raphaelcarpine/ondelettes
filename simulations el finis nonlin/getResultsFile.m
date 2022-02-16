function filePath = getResultsFile(fileNb, subFolder)

if exist('C:\Users\carpine\Documents\projets\simulations elements finis non lin\data', 'dir')
    dataFolder = 'C:\Users\carpine\Documents\projets\simulations elements finis non lin\data';
elseif exist('C:\Users\raphael\Documents\resultats simul diff finies', 'dir')
    dataFolder = 'C:\Users\raphael\Documents\resultats simul diff finies';
else
    error(' ');
end

if nargin > 1
    dataFolder = fullfile(dataFolder, subFolder);
end

% last simulation
if nargin == 0 || fileNb <= 0
    if nargin == 0
        fileNb = 0;
    end
    listing = dir(dataFolder);
    listingNames = {listing.name};
    nbSimul = 0;
    for kname = 1:length(listingNames)
        name1 = strsplit(listingNames{kname}, '_');
        name1 = name1{1};
        if length(name1) > 5 && strcmp(name1(1:5), 'simul')
            nbSimul = max(nbSimul, str2double(name1(6:end)));
        end
    end
    fileNb = nbSimul + fileNb;
end

% data download
listing = dir(dataFolder);
listingNames = {listing.name};
simul_name = '';
for kname = 1:length(listingNames)
    name1 = strsplit(listingNames{kname}, '_');
    name1 = name1{1};
    if length(name1) > 5 && strcmp(name1(1:5), 'simul') && str2double(name1(6:end)) == fileNb
        simul_name = listingNames{kname};
        break
    end
end
if isempty(simul_name)
    error('file not found');
end

filePath = fullfile(dataFolder, simul_name);

end

