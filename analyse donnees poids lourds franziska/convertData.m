fileFolder = 'C:\Users\carpine\Documents\projets\donnees poids lourds franziska';
fileName = '2015_14.csv';

T = readtable(fullfile(fileFolder, fileName), 'DatetimeType', 'text');
