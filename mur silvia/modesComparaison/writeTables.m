load('mur silvia\modesFourier\ModesLMS.mat');
load('mur silvia\modesOndelette\allModalQuantities.mat');
load('mur silvia\modesOndelette2\ModesWvlt2.mat');

allFreqs = AllModalQuantities.freqs;
allShapes = AllModalQuantities.shapes;
allDamps = AllModalQuantities.damps;

%%
createDoc = true;
showStdDevInp = false;

%%
P = [0, 6, 7];

nbModes = [3, 4, 3];

tableMatFreq = 0;
tableMatDamp = 0;
tableMatMAC = 0;
tableMatInp = 0;
lineMat = 1;
columnMatFreq = 1;
columnMatDamp = 1;
columnMatMAC = 1;
columnMatInp = 1;

%%

for indp = 1:3
    p = P(indp);
    
    for mode = 1:nbModes(indp)
        %% transients
        % nb transient
        nbTransients = length(allFreqs{indp}{mode});
        
        %% freqs
        % freq LMS
        tableMatFreq(lineMat, columnMatFreq) = ModesLMS(indp,mode).freq;
        columnMatFreq = columnMatFreq+1;
        
        % mean freqs cwt
        tableMatFreq(lineMat, columnMatFreq) = mean(allFreqs{indp}{mode});
        columnMatFreq = columnMatFreq+1;
        
        % freq CWT2
        if mode <= size(ModesWvlt2, 2) && ~isempty(ModesWvlt2(indp,mode).freq)
            tableMatFreq(lineMat, columnMatFreq) = ModesWvlt2(indp,mode).freq;
        else
            tableMatFreq(lineMat, columnMatFreq) = nan;
        end
        columnMatFreq = columnMatFreq+1;
        
        % nb transients
        tableMatFreq(lineMat, columnMatFreq) = nbTransients;
        columnMatFreq = columnMatFreq + 1;
        
        % std freqs cwt
        if nbTransients == 1
            tableMatFreq(lineMat, columnMatFreq) = nan;
        else
            tableMatFreq(lineMat, columnMatFreq) = std(allFreqs{indp}{mode});
        end
        columnMatFreq = columnMatFreq+1;
        
        %% damping
        % damping LMS
        tableMatDamp(lineMat, columnMatDamp) = 100*ModesLMS(indp,mode).damping;
        columnMatDamp = columnMatDamp+1;
        
        % mean dampings cwt
        tableMatDamp(lineMat, columnMatDamp) = 100*mean(allDamps{indp}{mode});
        columnMatDamp = columnMatDamp+1;
        
        % damping CWT2
        if mode <= size(ModesWvlt2, 2) && ~isempty(ModesWvlt2(indp,mode).damping)
            tableMatDamp(lineMat, columnMatDamp) = 100*ModesWvlt2(indp,mode).damping;
        else
            tableMatDamp(lineMat, columnMatDamp) = nan;
        end
        columnMatDamp = columnMatDamp+1;
        
        % nb transients
        tableMatDamp(lineMat, columnMatDamp) = nbTransients;
        columnMatDamp = columnMatDamp + 1;
        
        % std dampings cwtif nbTransients == 1
        if nbTransients == 1
            tableMatDamp(lineMat, columnMatDamp) = nan;
        else
            tableMatDamp(lineMat, columnMatDamp) = 100*std(allDamps{indp}{mode});
        end
        columnMatDamp = columnMatDamp+1;
        
        %% MAC
        shapeF = ModesLMS(indp,mode).shape;
        meanShape = transpose( mean( allShapes{indp}{mode}, 1));
        stdShape = std( real( allShapes{indp}{mode}), 0, 1) + 1i*std( imag( allShapes{indp}{mode}), 0, 1);
        if mode <= size(ModesWvlt2, 2) && ~isempty(ModesWvlt2(indp,mode).shape)
            shapeCWT2 = ModesWvlt2(indp,mode).shape;
        else
            shapeCWT2 = nan(9, 1);
        end
        
        % mac 12
        mac = abs(meanShape'*shapeF)^2 / ((meanShape'*meanShape) * (shapeF'*shapeF));
        tableMatMAC(lineMat, columnMatMAC) = 100*mac;
        columnMatMAC = columnMatMAC+1;
        
        % mac 13
        mac = abs(shapeCWT2'*shapeF)^2 / ((shapeCWT2'*shapeCWT2) * (shapeF'*shapeF));
        tableMatMAC(lineMat, columnMatMAC) = 100*mac;
        columnMatMAC = columnMatMAC+1;
        
        % mac 23
        mac = abs(shapeCWT2'*meanShape)^2 / ((shapeCWT2'*shapeCWT2) * (meanShape'*meanShape));
        tableMatMAC(lineMat, columnMatMAC) = 100*mac;
        columnMatMAC = columnMatMAC+1;
        
        %% Inp
        
        % I LMS
        tableMatInp(lineMat, columnMatInp) = 100*nonPropIndex(shapeF);
        columnMatInp = columnMatInp+1;
        
        % I cwt
        tableMatInp(lineMat, columnMatInp) = 100*nonPropIndex(meanShape);
        columnMatInp = columnMatInp+1;
        
        % I cwt2
        tableMatInp(lineMat, columnMatInp) = 100*nonPropIndex(shapeCWT2);
        columnMatInp = columnMatInp+1;
        
        if showStdDevInp
            % nb transients
            tableMatInp(lineMat, columnMatInp) = nbTransients;
            columnMatInp = columnMatInp + 1;
            
            % SD I
            % premiere methode
            allI = [];
            for ktransient = 1:size(allShapes{indp}{mode}, 1)
                allI = [allI, nonPropIndex( transpose( allShapes{indp}{mode}(ktransient, :)))];
            end
            stdI = std(allI);
            % deuxieme methode
            covI = cov(imag(allShapes{indp}{mode}));
            covR = cov(real(allShapes{indp}{mode}));
            covI = covI .* eye(9);
            covR = covR .* eye(9);
            stdI = 1/(nonPropIndex(meanShape) * norm(meanShape)^4) * sqrt(...
                norm(real(meanShape))^4 * imag(meanShape).' * covI * imag(meanShape) +...
                norm(imag(meanShape))^4 * real(meanShape).' * covR * real(meanShape) );
            
            if nbTransients == 1
                tableMatInp(lineMat, columnMatInp) = nan;
            else
                tableMatInp(lineMat, columnMatInp) = 100*stdI;
            end
            columnMatInp = columnMatInp+1;
        end
        
        %%
        columnMatFreq = 1;
        columnMatDamp = 1;
        columnMatMAC = 1;
        columnMatInp = 1;
        lineMat = lineMat+1;
        
    end
end







%% tex freq doc

tableFile = fopen('mur silvia\modesComparaison\templates\tableFreqTemplate.txt', 'r');
docStringFreq = fread(tableFile);
fclose(tableFile);

docStringFreq = char(transpose(docStringFreq));
docStringFreq = strrep(docStringFreq, 'Â', '');
docStringFreq = strrep(docStringFreq, '\%', 'µ');
docStringFreq = strrep(docStringFreq, '\', '§');
docStringFreq = strrep(docStringFreq, '£', '%');
docStringFreq = sprintf(docStringFreq, reshape(transpose(tableMatFreq), 1, []));
docStringFreq = strrep(docStringFreq, '§', '\');
docStringFreq = strrep(docStringFreq, 'µ', '\%');
docStringFreq = strrep(docStringFreq, 'NaN', '/');


%% tex damp doc

tableFile = fopen('mur silvia\modesComparaison\templates\tableDampTemplate.txt', 'r');
docStringDamp = fread(tableFile);
fclose(tableFile);

docStringDamp = char(transpose(docStringDamp));
docStringDamp = strrep(docStringDamp, 'Â', '');
docStringDamp = strrep(docStringDamp, '\%', 'µ');
docStringDamp = strrep(docStringDamp, '\', '§');
docStringDamp = strrep(docStringDamp, '£', '%');
docStringDamp = sprintf(docStringDamp, reshape(transpose(tableMatDamp), 1, []));
docStringDamp = strrep(docStringDamp, '§', '\');
docStringDamp = strrep(docStringDamp, 'µ', '\%');
docStringDamp = strrep(docStringDamp, 'NaN', '/');


%% tex MAC doc

tableFile = fopen('mur silvia\modesComparaison\templates\tableMACTemplate.txt', 'r');
docStringMAC = fread(tableFile);
fclose(tableFile);

docStringMAC = char(transpose(docStringMAC));
docStringMAC = strrep(docStringMAC, 'Â', '');
docStringMAC = strrep(docStringMAC, '\%', 'µ');
docStringMAC = strrep(docStringMAC, '\', '§');
docStringMAC = strrep(docStringMAC, '£', '%');
docStringMAC = sprintf(docStringMAC, reshape(transpose(tableMatMAC), 1, []));
docStringMAC = strrep(docStringMAC, '§', '\');
docStringMAC = strrep(docStringMAC, 'µ', '\%');
docStringMAC = strrep(docStringMAC, 'NaN', '/');


%% tex damp doc

if showStdDevInp
    tableFile = fopen('mur silvia\modesComparaison\templates\tableInpTemplate.txt', 'r');
else
    tableFile = fopen('mur silvia\modesComparaison\templates\tableInpTemplate2.txt', 'r');
end
docStringInp = fread(tableFile);
fclose(tableFile);

docStringInp = char(transpose(docStringInp));
docStringInp = strrep(docStringInp, 'Â', '');
docStringInp = strrep(docStringInp, '\%', 'µ');
docStringInp = strrep(docStringInp, '\', '§');
docStringInp = strrep(docStringInp, '£', '%');
docStringInp = sprintf(docStringInp, reshape(transpose(tableMatInp), 1, []));
docStringInp = strrep(docStringInp, '§', '\');
docStringInp = strrep(docStringInp, 'µ', '\%');
docStringInp = strrep(docStringInp, 'NaN', '/');



%% creation documents

% freq
tableFile = fopen('mur silvia\modesComparaison\tableFreq.tex', 'w');
fwrite(tableFile, docStringFreq);
fclose(tableFile);
% damp
tableFile = fopen('mur silvia\modesComparaison\tableDamp.tex', 'w');
fwrite(tableFile, docStringDamp);
fclose(tableFile);
% MAC
tableFile = fopen('mur silvia\modesComparaison\tableMAC.tex', 'w');
fwrite(tableFile, docStringMAC);
fclose(tableFile);
% Inp
tableFile = fopen('mur silvia\modesComparaison\tableInp.tex', 'w');
fwrite(tableFile, docStringInp);
fclose(tableFile);


%% creation document pour test
if createDoc
    headerString = ['\documentclass[11pt]{article}', newline, newline, newline, '\usepackage[margin=0.5in]{geometry}', newline, '\usepackage{multirow}', newline, newline, newline, '\begin{document}', newline, newline];
    endString = ['', newline, newline, '\end{document}'];
    headerTable = ['\begin{table}', newline, '\centering', newline];
    endTable = ['', newline, '\end{table}', newline, newline];
    
    tableDocFile = fopen('mur silvia\modesComparaison\document test 2\tableDoc.tex', 'w');
    fwrite(tableDocFile, [headerString,...
        headerTable, docStringFreq, endTable,...
        headerTable, docStringDamp, endTable,...
        headerTable, docStringMAC, endTable,...
        headerTable, docStringInp, endTable,...
        endString]);
    fclose(tableDocFile);
end
