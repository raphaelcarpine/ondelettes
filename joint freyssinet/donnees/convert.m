fileFolder = 'C:\Users\carpine\Documents\projets\donnees charles\';
saveFolder = 'joint freyssinet/donnees/';
filePath = [fileFolder, 'ZOOM0002_LR.WAV'];

%% enregistrement des donnees

% lecture des donnes
[X, Fs] = audioread(filePath);
X = transpose(X);

Nhalf = floor(size(X, 2)/2);
X11 = X(1, 1:Nhalf);
X12 = X(1, Nhalf+1:end);
X21 = X(2, 1:Nhalf);
X22 = X(2, Nhalf+1:end);

% enregistrement
save([saveFolder, 'data11'], 'X11', 'Fs');
save([saveFolder, 'data12'], 'X12', 'Fs');
save([saveFolder, 'data21'], 'X21', 'Fs');
save([saveFolder, 'data22'], 'X22', 'Fs');