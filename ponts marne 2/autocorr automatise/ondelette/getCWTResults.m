function [filesNames, FreqsCell, DampsCell, ShapesCell] = getCWTResults()

%% modes à ignorer (modes proches etc)

modIgnore = {
    [10],
    [1 2 10 11 12],
    [1 2 9 10 11],
    [1 2],
    [1 2],
    [1 2],
    [1 2],
    [1 2],
    [1 2 10]};


%% chargement data

dataPath = 'ponts marne 2\autocorr automatise\ondelette\save';

filesData = dir(fullfile(dataPath, '*.mat'));
filesPaths = {filesData.folder};
filesNames = {filesData.name};
for kfile = 1:length(filesNames)
    filesPaths{kfile} = fullfile(filesPaths{kfile}, filesNames{kfile});
    [~, filesNames{kfile}] = fileparts(filesNames{kfile});
end

ShapesCell = cell(1, length(filesNames));
FreqsCell = cell(1, length(filesNames));
DampsCell = cell(1, length(filesNames));

for kfile = 1:length(filesNames)
    load(filesPaths{kfile});
    I = true(1, length(Freqs));
    I(modIgnore{kfile}) = false;
    ShapesCell{kfile} = Shapes(:, I);
    FreqsCell{kfile} = Freqs(I);
    DampsCell{kfile} = Damps(I);
end