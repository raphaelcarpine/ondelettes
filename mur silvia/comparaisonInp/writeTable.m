load('mur silvia\modesFourier\ModesLMS.mat');
load('mur silvia\modesOndelette\allModalQuantities.mat');
load('mur silvia\modesOndelette2\ModesWvlt2.mat');

allFreqs = AllModalQuantities.freqs;
allShapes = AllModalQuantities.shapes;
allDamps = AllModalQuantities.damps;

FreqsAdhikari = 1/2/pi * [32.5838, 33.5570];
ZetaAdhikari = 100 * [1/2/11.8178, 1/2/6.4006];
InpAdhikari = [49.26, 51.68];
InpTildeAdhikari = [41.19, 42.75];

FreqsErlicher = [26.88, 35.26, 59.53, 65.58, 75.92, 107.08, 128.02, 129.55, 158.57, 179.46];
ZetaErlicher = [2.7, 3.0, 2.8, 3.1, 3.05, 2.95, 1.80, 2.9, 2.2, 2.35];
InpTildeErlicher = 100 * [0.072, 0.309, 0.049, 0.042, 0.077, 0.175, 0.471, 0.337, 0.1465, 0.107];

%%
createDoc = true;
CWT1only = false;
excludeCWT2 = true;

%%
P = [0, 6, 7];

nbModes = [3, 4, 3];

%% calcul moyennes

FreqsWall = nan(3, 4);
ZetaWall = nan(3, 4);
InpWall = nan(3, 4);
for kp = 1:3
    for kmode = 1:nbModes(kp)
        if CWT1only
            FreqsWall(kp, kmode) = mean(allFreqs{kp}{kmode});
            ZetaWall(kp, kmode) = mean(allDamps{kp}{kmode});
            InpWall(kp, kmode) = nonPropIndex( transpose( mean( allShapes{kp}{kmode}, 1)));
        else
            FreqsArray = mean(allFreqs{kp}{kmode});
            ZetaArray = mean(allDamps{kp}{kmode});
            InpArray = nonPropIndex( transpose( mean( allShapes{kp}{kmode}, 1)));
            
            FreqsArray(end+1) = ModesLMS(kp, kmode).freq;
            ZetaArray(end+1) = ModesLMS(kp, kmode).damping;
            if ~(kp == 3 && kmode == 3) % valeur écartée pour LSCF
                InpArray(end+1) = nonPropIndex(ModesLMS(kp, kmode).shape);
            end
            
            if ~excludeCWT2 && kmode <= 3 && ~isempty(ModesWvlt2(kp, kmode).freq)
                FreqsArray(end+1) = ModesWvlt2(kp, kmode).freq;
                ZetaArray(end+1) = ModesWvlt2(kp, kmode).damping;
                InpArray(end+1) = nonPropIndex(ModesWvlt2(kp, kmode).shape);
            end
            
            FreqsWall(kp, kmode) = mean(FreqsArray);
            ZetaWall(kp, kmode) = mean(ZetaArray);
            InpWall(kp, kmode) = mean(InpArray);
        end
    end
end



%% vecteur valeurs

dataVect = [];

% wall
for kp = 1:3
    for kmode = 1:nbModes(kp)
        if kmode > 1
            eta = 2*abs(FreqsWall(kp, kmode) - FreqsWall(kp, kmode-1)) / (FreqsWall(kp, kmode) + FreqsWall(kp, kmode-1));
            dataVect(end+1) = 100 * eta;
        end
        dataVect(end+1) = FreqsWall(kp, kmode);
        dataVect(end+1) = 100 * ZetaWall(kp, kmode);
        dataVect(end+1) = nan;
        dataVect(end+1) = 100 * InpWall(kp, kmode);
    end
end

% adhikari
for kmode = 1:length(FreqsAdhikari)
    if kmode > 1
        eta = 2*abs(FreqsAdhikari(kmode) - FreqsAdhikari(kmode-1)) / (FreqsAdhikari(kmode) + FreqsAdhikari(kmode-1));
        dataVect(end+1) = 100 * eta;
    end
    dataVect(end+1) = FreqsAdhikari(kmode);
    dataVect(end+1) = ZetaAdhikari(kmode);
    dataVect(end+1) = InpAdhikari(kmode);
    dataVect(end+1) = InpTildeAdhikari(kmode);
end

% erlicher
for kmode = 1:length(FreqsErlicher)
    if kmode > 1
        eta = 2*abs(FreqsErlicher(kmode) - FreqsErlicher(kmode-1)) / (FreqsErlicher(kmode) + FreqsErlicher(kmode-1));
        dataVect(end+1) = 100 * eta;
    end
    dataVect(end+1) = FreqsErlicher(kmode);
    dataVect(end+1) = ZetaErlicher(kmode);
        dataVect(end+1) = nan;
    dataVect(end+1) = InpTildeErlicher(kmode);
end

%% tex doc

tableFile = fopen('mur silvia\comparaisonInp\template\tableTemplate.txt', 'r');
docString = fread(tableFile);
fclose(tableFile);

docString = char(transpose(docString));
docString = strrep(docString, 'Â', '');
docString = strrep(docString, '\%', 'µ');
docString = strrep(docString, '%', '¤');
docString = strrep(docString, '\', '§');
docString = strrep(docString, '£', '%');
docString = sprintf(docString, dataVect);
docString = strrep(docString, '§', '\');
docString = strrep(docString, '¤', '%');
docString = strrep(docString, 'µ', '\%');
docString = strrep(docString, 'NaN', '/');




%% creation document

tableFile = fopen('mur silvia\comparaisonInp\tableComparison.tex', 'w');
fwrite(tableFile, docString);
fclose(tableFile);


%% creation document pour test
if createDoc
    headerString = ['\documentclass[11pt]{article}', newline, newline, newline, '\usepackage[margin=0.5in]{geometry}', newline,...
        '\usepackage{multirow}', newline, newline, '\renewcommand{\tabcolsep}{10 pt}', newline,...
        '\renewcommand{\arraystretch}{.6}', newline, '\usepackage{natbib}',...
        newline, newline, newline, '\begin{document}', newline, newline];
    endString = ['', newline, newline, '\end{document}'];
    headerTable = ['\begin{table}', newline, '\centering', newline];
    endTable = ['', newline, '\end{table}', newline, newline];
    
    tableDocFile = fopen('mur silvia\comparaisonInp\document test\tableDoc.tex', 'w');
    fwrite(tableDocFile, [headerString, headerTable, docString, endTable, endString]);
    fclose(tableDocFile);
end
