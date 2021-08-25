function [X, t, labels] = getDataZ24(measure, config)
%GETDATA Summary of this function goes here
%   mes: measure
%   conf: configuration

if measure < 1 || measure > 17
    error('1 <= mes <= 17');
end
if config < 1 || config > 9
    error('1 <= conf <= 9');
end

%% data

dt = 0.01;
conversionCoeff = 5; % V/g
g = 9.81;

%% file path

dataFolder = 'C:\Users\carpine\Documents\projets\Z24\data';

if measure < 9
    subFolder = 'pdt_01_08';
else
    subFolder = 'pdt_09_17';
end

subsubFolder = sprintf('%02u', measure);

subsubsubFolder = 'avt'; % ambiant vibration test

fileName = sprintf('%02usetup%02u.mat', [measure, config]);

%% loading & formattting

filePath = fullfile(dataFolder, subFolder, subsubFolder, subsubsubFolder, fileName);

load(filePath);

X = transpose(data);
t = (0:size(X, 2)-1) * dt;

labels = {};
for k = 1:size(labelshulp, 1)
    label = labelshulp(k, :);
    while ~isempty(label) && label(end) == ' '
        label = label(1:end-1);
    end
    labels{end+1} = label;
end

%% conversion

X = (g/conversionCoeff)*X;

end

