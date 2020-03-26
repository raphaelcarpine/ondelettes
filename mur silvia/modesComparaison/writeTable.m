load('mur silvia\modesFourier\ModesLMS.mat');


load('mur silvia\modesOndelette\allModalQuantities.mat');

allFreqs = AllModalQuantities.freqs;
allShapes = AllModalQuantities.shapes;
allDamps = AllModalQuantities.damps;

%%
createDoc = true;

%%
P = [0, 6, 7];

tableMat = 0;
lineMat = 1;
columnMat = 1;

%%

for indp = 1:3
    p = P(indp);
    
    for mode = 1:nbModes(indp)
        
        %% transients
        % nb transient
        nbTransients = length(allFreqs{indp}{mode});
        tableMat(lineMat, columnMat) = nbTransients;
        columnMat = columnMat+1;
        
        %% freqs
        % freq LMS
        tableMat(lineMat, columnMat) = ModesLMS(indp,mode).freq;
        columnMat = columnMat+1;
        
        % mean freqs cwt
        tableMat(lineMat, columnMat) = mean(allFreqs{indp}{mode});
        columnMat = columnMat+1;
        
        % std freqs cwt
        if nbTransients == 1
            tableMat(lineMat, columnMat) = nan;
        else
            tableMat(lineMat, columnMat) = std(allFreqs{indp}{mode});
        end
        columnMat = columnMat+1;
        
        %% damping
        % damping LMS
        tableMat(lineMat, columnMat) = 100*ModesLMS(indp,mode).damping;
        columnMat = columnMat+1;
        
        % mean dampings cwt
        tableMat(lineMat, columnMat) = 100*mean(allDamps{indp}{mode});
        columnMat = columnMat+1;
        
        % std dampings cwtif nbTransients == 1
        if nbTransients == 1
            tableMat(lineMat, columnMat) = nan;
        else
            tableMat(lineMat, columnMat) = 100*std(allDamps{indp}{mode});
        end
        columnMat = columnMat+1;
        
        %% shape
        meanShape = transpose( mean( allShapes{indp}{mode}, 1));
        stdShape = std( real( allShapes{indp}{mode}), 0, 1) + 1i*std( imag( allShapes{indp}{mode}), 0, 1);
        shapeF = ModesLMS(indp,mode).shape;
        
        % mac
        mac = abs(meanShape'*shapeF)^2 / ((meanShape'*meanShape) * (shapeF'*shapeF));
        tableMat(lineMat, columnMat) = 100*mac;
        columnMat = columnMat+1;
        
        % I LMS
        tableMat(lineMat, columnMat) = 100*nonPropIndex(shapeF);
        columnMat = columnMat+1;
        
        % I cwt
        tableMat(lineMat, columnMat) = 100*nonPropIndex(meanShape);
        columnMat = columnMat+1;
        
        % SD I
        allI = [];
        for ktransient = 1:size(allShapes{indp}{mode}, 1)
            allI = [allI, nonPropIndex( transpose( allShapes{indp}{mode}(ktransient, :)))];
        end
        if nbTransients == 1
            tableMat(lineMat, columnMat) = nan;
        else
            tableMat(lineMat, columnMat) = 100*std( allI);
        end
        columnMat = columnMat+1;
        
        %%
        columnMat = 1;
        lineMat = lineMat+1;
        
    end
end









%% tex doc

docString = '';
docString = [docString, '\begin{tabular}{ccc|c|c|c|c|c|c|c|c|c|c|} ', newline];
docString = [docString, '\cline{4-13} &  &  & \multicolumn{3}{c|}{Frequencies} & \multicolumn{3}{c|}{Damping coefficients}  & \multicolumn{4}{c|}{Modal shapes} \\ \hline ', newline];
docString = [docString, '\multicolumn{1}{|c|}{Step} & ', newline];
docString = [docString, '\multicolumn{1}{c|}{\begin{tabular}[c]{@{}c@{}} Mode\\ number \end{tabular}} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}Number of\\ transients\\ (CWT) \end{tabular} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}Frequency\\ (LSCF)\\ Hz \end{tabular} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}Frequency\\ (CWT)\\ Hz\end{tabular} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}Std Dev\\ (CWT)\\ Hz\end{tabular} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}Damping\\ (LSCF)\\ \% \end{tabular} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}Damping\\ (CWT)\\ \% \end{tabular} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}Std Dev\\ (CWT)\\ \% \end{tabular} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}MAC\\ (LMS\\ $\times$ CWT)\\ \% \end{tabular} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}$\tilde I_{np}$\\ (LSCF)\\ \% \end{tabular} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}$\tilde I_{np}$\\ (CWT)\\ \% \end{tabular} & ', newline];
docString = [docString, '\begin{tabular}[c]{@{}c@{}}Std Dev\\ (CWT)\\ \% \end{tabular} ', newline];
docString = [docString, ' \\ \hline ', newline, newline];


lineMat = 1;
for indp = 1:3
    p = P(indp);
    
    for mode = 1:nbModes(indp)
        if mode == 1
            docString = [docString, '\multicolumn{1}{|c|}{\multirow{3}{*}{$P_', num2str(p), '$}} & \multicolumn{1}{c|}{1} ', newline];
        else
            docString = [docString, '\multicolumn{1}{|c|}{} & \multicolumn{1}{c|}{', num2str(mode), '} ', newline];
        end
        
        for columnMat = 1:size(tableMat, 2)
            if isnan( tableMat(lineMat, columnMat))
                docString = [docString, ' & /'];
            elseif columnMat == 1
                docString = [docString, ' & ', num2str(tableMat(lineMat, columnMat), '%u')];
            else
                docString = [docString, ' & ', num2str(tableMat(lineMat, columnMat), '%.2f')];
            end
        end
        
        if mode ~= nbModes(indp)
            docString = [docString, ' \\ \cline{2-13} ', newline];
        else
            docString = [docString, ' \\ \hline ', newline, newline];
        end
        
        lineMat = lineMat+1;
    end
end

docString = [docString, '\end{tabular}'];


%% creation document

tableFile = fopen('mur silvia\modesComparaison\table.txt', 'w');
fwrite(tableFile, docString);
fclose(tableFile);


%% creation document pour test
if createDoc
    headerString = ['\documentclass[11pt, landscape]{article}', newline, newline, '\usepackage[margin=0.1in]{geometry}', newline, newline, '\usepackage{multirow}', newline, '\usepackage{graphicx}', newline, newline, '\begin{document}', newline, newline, '\begin{small}', newline, newline, '\begin{table}', newline];
    endString = ['', newline, '\end{table}', newline, newline, '\end{small}', newline, newline, '\end{document}'];
    
    tableDocFile = fopen('mur silvia\modesComparaison\document test\tableDoc.tex', 'w');
    fwrite(tableDocFile, [headerString, docString, endString]);
    fclose(tableDocFile);
end
