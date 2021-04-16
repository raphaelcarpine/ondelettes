fileFolder = 'C:\Users\carpine\Documents\projets\donnees charles';
saveFolder = fileFolder;

%% enregistrement des donnees

files = dir(fullfile(fileFolder, '*.wav'));

for k=1:length(files)
    file = files(k);
    
    % lecture des donnes
    [X, Fs] = audioread(fullfile(fileFolder, file.name));
    X = transpose(X);
    
    % enregistrement
    [~, nameFile] = fileparts(file.name);
    save(fullfile(saveFolder, nameFile), 'X', 'Fs');
end