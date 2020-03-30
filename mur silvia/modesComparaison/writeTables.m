load('mur silvia\modesFourier\ModesLMS.mat');
load('mur silvia\modesOndelette\allModalQuantities.mat');
load('mur silvia\modesOndelette2\ModesWvlt2.mat');

allFreqs = AllModalQuantities.freqs;
allShapes = AllModalQuantities.shapes;
allDamps = AllModalQuantities.damps;

%%
createDoc = true;

%%
P = [0, 6, 7];

nbModes = [3, 4, 3];

tableMatFreq = 0;
tableMatDamp = 0;
tableMatShape = 0;
lineMat = 1;
columnMatFreq = 1;
columnMatDamp = 1;
columnMatShape = 1;

%%

for indp = 1:3
    p = P(indp);
    
    for mode = 1:nbModes(indp)
        
        %% transients
        % nb transient
        nbTransients = length(allFreqs{indp}{mode});
        tableMatFreq(lineMat, columnMatFreq) = nbTransients;
        tableMatDamp(lineMat, columnMatDamp) = nbTransients;
        tableMatShape(lineMat, columnMatShape) = nbTransients;
        columnMatFreq = columnMatFreq + 1;
        columnMatDamp = columnMatDamp + 1;
        columnMatShape = columnMatShape + 1;
        
        %% freqs
        % freq LMS
        tableMatFreq(lineMat, columnMatFreq) = ModesLMS(indp,mode).freq;
        columnMatFreq = columnMatFreq+1;
        
        % mean freqs cwt
        tableMatFreq(lineMat, columnMatFreq) = mean(allFreqs{indp}{mode});
        columnMatFreq = columnMatFreq+1;
        
        % std freqs cwt
        if nbTransients == 1
            tableMatFreq(lineMat, columnMatFreq) = nan;
        else
            tableMatFreq(lineMat, columnMatFreq) = std(allFreqs{indp}{mode});
        end
        columnMatFreq = columnMatFreq+1;
        
        % freq CWT2
        if mode <= size(ModesWvlt2, 2) && ~isempty(ModesWvlt2(indp,mode).freq)
            tableMatFreq(lineMat, columnMatFreq) = ModesWvlt2(indp,mode).freq;
        else
            tableMatFreq(lineMat, columnMatFreq) = nan;
        end
        columnMatFreq = columnMatFreq+1;
        
        %% damping
        % damping LMS
        tableMatDamp(lineMat, columnMatDamp) = 100*ModesLMS(indp,mode).damping;
        columnMatDamp = columnMatDamp+1;
        
        % mean dampings cwt
        tableMatDamp(lineMat, columnMatDamp) = 100*mean(allDamps{indp}{mode});
        columnMatDamp = columnMatDamp+1;
        
        % std dampings cwtif nbTransients == 1
        if nbTransients == 1
            tableMatDamp(lineMat, columnMatDamp) = nan;
        else
            tableMatDamp(lineMat, columnMatDamp) = 100*std(allDamps{indp}{mode});
        end
        columnMatDamp = columnMatDamp+1;
        
        % damping CWT2
        if mode <= size(ModesWvlt2, 2) && ~isempty(ModesWvlt2(indp,mode).damping)
            tableMatDamp(lineMat, columnMatDamp) = 100*ModesWvlt2(indp,mode).damping;
        else
            tableMatDamp(lineMat, columnMatDamp) = nan;
        end
        columnMatDamp = columnMatDamp+1;
        
        %% shape
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
        tableMatShape(lineMat, columnMatShape) = 100*mac;
        columnMatShape = columnMatShape+1;
        
        % mac 13
        mac = abs(shapeCWT2'*shapeF)^2 / ((shapeCWT2'*shapeCWT2) * (shapeF'*shapeF));
        tableMatShape(lineMat, columnMatShape) = 100*mac;
        columnMatShape = columnMatShape+1;
        
        % mac 23
        mac = abs(shapeCWT2'*meanShape)^2 / ((shapeCWT2'*shapeCWT2) * (meanShape'*meanShape));
        tableMatShape(lineMat, columnMatShape) = 100*mac;
        columnMatShape = columnMatShape+1;
        
        % I LMS
        tableMatShape(lineMat, columnMatShape) = 100*nonPropIndex(shapeF);
        columnMatShape = columnMatShape+1;
        
        % I cwt
        tableMatShape(lineMat, columnMatShape) = 100*nonPropIndex(meanShape);
        columnMatShape = columnMatShape+1;
        
        % SD I
        allI = [];
        for ktransient = 1:size(allShapes{indp}{mode}, 1)
            allI = [allI, nonPropIndex( transpose( allShapes{indp}{mode}(ktransient, :)))];
        end
        if nbTransients == 1
            tableMatShape(lineMat, columnMatShape) = nan;
        else
            tableMatShape(lineMat, columnMatShape) = 100*std( allI);
        end
        columnMatShape = columnMatShape+1;
        
        % I cwt2
        tableMatShape(lineMat, columnMatShape) = 100*nonPropIndex(shapeCWT2);
        columnMatShape = columnMatShape+1;
        
        %%
        columnMatFreq = 1;
        columnMatDamp = 1;
        columnMatShape = 1;
        lineMat = lineMat+1;
        
    end
end







%% tex freq doc

docStringFreq = ['\begin{tabular}{cc', repmat('c|', 1, size(tableMatFreq, 2)), '} ', newline];
docStringFreq = [docStringFreq, '\cline{4-', num2str(size(tableMatFreq, 2)+2), '} &  &  & \multicolumn{', num2str(size(tableMatFreq, 2)-1), '}{c|}{Frequencies} \\ \hline ', newline];
docStringFreq = [docStringFreq, '\multicolumn{1}{|c|}{Step} & ', newline];
docStringFreq = [docStringFreq, '\multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}} Mode\\ number \end{tabular}} & ', newline];
docStringFreq = [docStringFreq, '\begin{tabular}[c]{@{}c@{}}Number of\\ transients\\ (CWT) \end{tabular} & ', newline];
docStringFreq = [docStringFreq, '\begin{tabular}[c]{@{}c@{}}Frequency\\ (LSCF)\\ Hz \end{tabular} & ', newline];
docStringFreq = [docStringFreq, '\begin{tabular}[c]{@{}c@{}}Frequency\\ (CWT)\\ Hz \end{tabular} & ', newline];
docStringFreq = [docStringFreq, '\begin{tabular}[c]{@{}c@{}}Std Dev\\ (CWT)\\ Hz \end{tabular} & ', newline];
docStringFreq = [docStringFreq, '\begin{tabular}[c]{@{}c@{}}Frequency\\ (CWT2)\\ Hz\end{tabular}', newline];
docStringFreq = [docStringFreq, ' \\ \hline  \hline', newline, newline];


lineMat = 1;
for indp = 1:3
    p = P(indp);
    
    for mode = 1:nbModes(indp)
        if mode == 1
            docStringFreq = [docStringFreq, '\multicolumn{1}{|c|}{\multirow{3}{*}{$P_', num2str(p), '$}} & \multicolumn{1}{c|}{1} ', newline];
        else
            docStringFreq = [docStringFreq, '\multicolumn{1}{|c|}{} & \multicolumn{1}{c|}{', num2str(mode), '} ', newline];
        end
        
        for columnMat = 1:size(tableMatFreq, 2)
            if isnan( tableMatFreq(lineMat, columnMat))
                docStringFreq = [docStringFreq, ' & /'];
            elseif columnMat == 1
                docStringFreq = [docStringFreq, ' & ', num2str(tableMatFreq(lineMat, columnMat), '%u')];
            else
                docStringFreq = [docStringFreq, ' & ', num2str(tableMatFreq(lineMat, columnMat), '%.2f')];
            end
        end
        
        if mode ~= nbModes(indp)
            docStringFreq = [docStringFreq, ' \\ \cline{2-', num2str(size(tableMatFreq, 2)+2), '} ', newline];
        elseif indp ~= 3
            docStringFreq = [docStringFreq, ' \\ \hline \hline', newline, newline];
        else
            docStringFreq = [docStringFreq, ' \\ \hline', newline, newline];
        end
        
        lineMat = lineMat+1;
    end
end

docStringFreq = [docStringFreq, '\end{tabular}'];

%% tex damp doc

docStringDamp = ['\begin{tabular}{cc', repmat('c|', 1, size(tableMatDamp, 2)), '} ', newline];
docStringDamp = [docStringDamp, '\cline{4-', num2str(size(tableMatDamp, 2)+2), '} &  &  & \multicolumn{', num2str(size(tableMatDamp, 2)-1), '}{c|}{Dampings} \\ \hline ', newline];
docStringDamp = [docStringDamp, '\multicolumn{1}{|c|}{Step} & ', newline];
docStringDamp = [docStringDamp, '\multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}} Mode\\ number \end{tabular}} & ', newline];
docStringDamp = [docStringDamp, '\begin{tabular}[c]{@{}c@{}}Number of\\ transients\\ (CWT) \end{tabular} & ', newline];
docStringDamp = [docStringDamp, '\begin{tabular}[c]{@{}c@{}}Damping\\ (LSCF)\\ \% \end{tabular} & ', newline];
docStringDamp = [docStringDamp, '\begin{tabular}[c]{@{}c@{}}Damping\\ (CWT)\\ \% \end{tabular} & ', newline];
docStringDamp = [docStringDamp, '\begin{tabular}[c]{@{}c@{}}Std Dev\\ (CWT)\\ \% \end{tabular} & ', newline];
docStringDamp = [docStringDamp, '\begin{tabular}[c]{@{}c@{}}Damping\\ (CWT2)\\ \% \end{tabular}', newline];
docStringDamp = [docStringDamp, ' \\ \hline \hline', newline, newline];


lineMat = 1;
for indp = 1:3
    p = P(indp);
    
    for mode = 1:nbModes(indp)
        if mode == 1
            docStringDamp = [docStringDamp, '\multicolumn{1}{|c|}{\multirow{3}{*}{$P_', num2str(p), '$}} & \multicolumn{1}{c|}{1} ', newline];
        else
            docStringDamp = [docStringDamp, '\multicolumn{1}{|c|}{} & \multicolumn{1}{c|}{', num2str(mode), '} ', newline];
        end
        
        for columnMat = 1:size(tableMatDamp, 2)
            if isnan( tableMatDamp(lineMat, columnMat))
                docStringDamp = [docStringDamp, ' & /'];
            elseif columnMat == 1
                docStringDamp = [docStringDamp, ' & ', num2str(tableMatDamp(lineMat, columnMat), '%u')];
            else
                docStringDamp = [docStringDamp, ' & ', num2str(tableMatDamp(lineMat, columnMat), '%.2f')];
            end
        end
        
        if mode ~= nbModes(indp)
            docStringDamp = [docStringDamp, ' \\ \cline{2-', num2str(size(tableMatDamp, 2)+2), '} ', newline];
        elseif indp ~= 3
            docStringDamp = [docStringDamp, ' \\ \hline \hline', newline, newline];
        else
            docStringDamp = [docStringDamp, ' \\ \hline', newline, newline];
        end
        
        lineMat = lineMat+1;
    end
end

docStringDamp = [docStringDamp, '\end{tabular}'];

%% tex shape doc

docStringShape = ['\begin{tabular}{cc', repmat('c|', 1, size(tableMatShape, 2)), '} ', newline];
docStringShape = [docStringShape, '\cline{4-', num2str(size(tableMatShape, 2)+2), '} &  &  & \multicolumn{', num2str(size(tableMatShape, 2)-1), '}{c|}{Modal Shapes} \\ \hline ', newline];
docStringShape = [docStringShape, '\multicolumn{1}{|c|}{Step} & ', newline];
docStringShape = [docStringShape, '\multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}} Mode\\ number \end{tabular}} & ', newline];
docStringShape = [docStringShape, '\begin{tabular}[c]{@{}c@{}}Number of\\ transients\\ (CWT) \end{tabular} & ', newline];
docStringShape = [docStringShape, '\begin{tabular}[c]{@{}c@{}}MAC\\ (CWT\\ $\times$ LSCF)\\ \% \end{tabular} & ', newline];
docStringShape = [docStringShape, '\begin{tabular}[c]{@{}c@{}}MAC\\ (CWT2\\ $\times$ LSCF)\\ \% \end{tabular} & ', newline];
docStringShape = [docStringShape, '\begin{tabular}[c]{@{}c@{}}MAC\\ (CWT2\\ $\times$ CWT)\\ \% \end{tabular} & ', newline];
docStringShape = [docStringShape, '\begin{tabular}[c]{@{}c@{}}$\tilde I_{np}$\\ (LSCF)\\ \% \end{tabular} & ', newline];
docStringShape = [docStringShape, '\begin{tabular}[c]{@{}c@{}}$\tilde I_{np}$\\ (CWT)\\ \% \end{tabular} & ', newline];
docStringShape = [docStringShape, '\begin{tabular}[c]{@{}c@{}}Std Dev\\ (CWT)\\ \% \end{tabular} & ', newline];
docStringShape = [docStringShape, '\begin{tabular}[c]{@{}c@{}}$\tilde I_{np}$\\ (CWT2)\\ \% \end{tabular}', newline];
docStringShape = [docStringShape, ' \\ \hline \hline', newline, newline];


lineMat = 1;
for indp = 1:3
    p = P(indp);
    
    for mode = 1:nbModes(indp)
        if mode == 1
            docStringShape = [docStringShape, '\multicolumn{1}{|c|}{\multirow{3}{*}{$P_', num2str(p), '$}} & \multicolumn{1}{c|}{1} ', newline];
        else
            docStringShape = [docStringShape, '\multicolumn{1}{|c|}{} & \multicolumn{1}{c|}{', num2str(mode), '} ', newline];
        end
        
        for columnMat = 1:size(tableMatShape, 2)
            if isnan( tableMatShape(lineMat, columnMat))
                docStringShape = [docStringShape, ' & /'];
            elseif columnMat == 1
                docStringShape = [docStringShape, ' & ', num2str(tableMatShape(lineMat, columnMat), '%u')];
            else
                docStringShape = [docStringShape, ' & ', num2str(tableMatShape(lineMat, columnMat), '%.2f')];
            end
        end
        
        if mode ~= nbModes(indp)
            docStringShape = [docStringShape, ' \\ \cline{2-', num2str(size(tableMatShape, 2)+2), '} ', newline];
        elseif indp ~= 3
            docStringShape = [docStringShape, ' \\ \hline \hline', newline, newline];
        else
            docStringShape = [docStringShape, ' \\ \hline', newline, newline];
        end
        
        lineMat = lineMat+1;
    end
end

docStringShape = [docStringShape, '\end{tabular}'];


%% creation documents

% freq
tableFile = fopen('mur silvia\modesComparaison\tableFreq.tex', 'w');
fwrite(tableFile, docStringFreq);
fclose(tableFile);
% damp
tableFile = fopen('mur silvia\modesComparaison\tableDamp.tex', 'w');
fwrite(tableFile, docStringDamp);
fclose(tableFile);
% shape
tableFile = fopen('mur silvia\modesComparaison\tableShape.tex', 'w');
fwrite(tableFile, docStringShape);
fclose(tableFile);


%% creation document pour test
if createDoc
    headerString = ['\documentclass[11pt]{article}', newline, newline, newline, '\usepackage[margin=0.5in]{geometry}', newline, '\usepackage{multirow}', newline, newline, newline, '\begin{document}', newline, newline];
    endString = ['', newline, newline, '\end{document}'];
    headerTable = ['\begin{table}', newline];
    endTable = ['', newline, '\end{table}', newline, newline];
    
    tableDocFile = fopen('mur silvia\modesComparaison\document test 2\tableDoc.tex', 'w');
    fwrite(tableDocFile, [headerString,...
        headerTable, docStringFreq, endTable,...
        headerTable, docStringDamp, endTable,...
        headerTable, docStringShape, endTable,...
        endString]);
    fclose(tableDocFile);
end
