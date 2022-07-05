function [filePath, fileName, filePathRadio, fileNameRadio, Nlag] = choixData(n)

if nargin == 0
    n = 0;
end

%%
dataFolder = 'C:\Users\carpine\Documents\projets\ponts marne\reprise operations 2022\donnees';
files = dir(fullfile(dataFolder, 'acceleros', '*.mat'));
files = {files.name};
filesNames = {};

%% removes "niveau" files & "error" files

k = 1;
while k <=length(files)
    if contains(files{k}, 'niveau')
        files = files([1:k-1, k+1:end]);
    elseif contains(files{k}, 'error')
        files = files([1:k-1, k+1:end]);
    else
        [~, filesNames{end+1}] = fileparts(files{k}); % removes extension
        k = k+1;
    end
end

%% sort by date

% filesDates = nan(1, length(filesNames));
% for k = 1:length(filesNames)
%     nameParts = strsplit(filesNames{k}, '_');
%     filesDates(k) = str2double(nameParts{2});
%     filesDates(k) = 100*mod(filesDates(k), 100) + floor(filesDates(k)/100); % month/day swap
% end
% [~, I] = sort(filesDates);
% 
% files = files(I);
% filesNames = filesNames(I);

%% choice

if n == 0
    % dialog box
    mlock
    persistent fileIndex;
    if isempty(fileIndex)
        fileIndex = 1;
    end
    
    fileIndex = listdlg('ListSize', [200 200], 'ListString', filesNames, 'SelectionMode', 'single', 'InitialValue', fileIndex); % 'InitialValue', FourierModeShapeMenu.UserData,...
    
else
    fileIndex = n;
end


%%
filePath = fullfile(dataFolder, 'acceleros', files{fileIndex});
fileName = filesNames{fileIndex};
filePathRadio = fullfile(dataFolder, 'radio', [filesNames{fileIndex}, '_resampled_filtered.mat']);
fileNameRadio = [filesNames{fileIndex}, '_resampled_filtered'];

Nlag = nan;
switch fileName
    case 'changis3006'
        Nlag = 162;
    case 'esbly0107_matin'
        Nlag = 416;
    case 'esbly0107_aprem'
        Nlag = 1750;
    otherwise
        error(' ');
end

end