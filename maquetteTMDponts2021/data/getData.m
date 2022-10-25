dataFolder = 'C:\Users\carpine\Documents\cours\semaine TMD ponts 2021\TP1\data';
% dataFile = 'SensorConnectDataTest16.csv';

%% position capteurs
% maquette
H = 212; % hauteur totale
h = 88; % hauteur de la poutre
l = 302.5; % demie longueur de la poutre

captNumber = [
    40196  40198  40199  40200  40201  40202  ];
captPos = [
    0      -l     -l/2   0      l/2    l      ; % x
    H      h      h      h      h      h      ]; % z

% captNumber = [
%     40200 40201  ];
% captNumber = 40201;


%% data

D = readtable(fullfile(dataFolder, dataFile), 'VariableNamingRule', 'preserve'); % conversion de 'name'.csv en tableau matlab 'X'
X = table2array(D(:, 2:end));
T = table2array(D(:, 1));
T = T - T(1);
T = seconds(T);
channelNames = D.Properties.VariableNames(2:end);
startDate = D.Time(1);
startDate.Year = startDate.Year + 2000;
startDate.TimeZone = 'UTC';
startDate.TimeZone = 'Europe/Paris';

X = X.';
T = T.';

%% data processing

[X, T] = removeRedundantData(X, T);
[X, T] = removeNanSignal(X, T);

[X, channelNames, captNumber2, captDir] = convertChNames(X, channelNames);

% if ~isequal(captNumber, captNumber2)
%     error('~isequal(captNumber, captNumber2)');
% end
