xlsxPath = 'C:\Users\carpine\Documents\projets\article pont de sens\tables\ResultsABharmonic.xlsx';
texPath = 'C:\Users\carpine\Documents\projets\article pont de sens\tables\ResultsABharmonic.tex';

formatSpec = {'', '%.2f', '%.1f', '%.2f', '%.2f', '%.2f', '%.1f', '%.2f', '%.1f', '%.2f', '%.1f'};
columnsSpec = '';

%% conversion xlsx

C = readcell(xlsxPath);

%% conversion en texte

for i = 1:size(C, 1)
    for j = 1:size(C, 2)
        if isnumeric(C{i, j})
            if isempty(formatSpec)
                C{i, j} = num2str(C{i, j});
            elseif length(formatSpec) == 1
                C{i, j} = sprintf(formatSpec{1}, C{i, j});
            elseif isempty(formatSpec{j})
                C{i, j} = num2str(C{i, j});
            else
                C{i, j} = sprintf(formatSpec{j}, C{i, j});
            end
        elseif ischar(C{i, j})
            % ok
        elseif ismissing(C{i, j})
            C{i, j} = '';
        end
    end
end

% egalisation longueurs char
for j = 1:size(C, 2)
    lmax = 0;
    for i = 1:size(C, 1)
        lmax = max(lmax, length(C{i, j}));
    end
    for i = 1:size(C, 1)
        C{i, j} = [C{i, j}, repmat(' ', 1, lmax - length(C{i, j}))];
    end
end


%% conversion en tex

s = '\begin{tabular}{';
if isempty(columnsSpec)
    s = [s, repmat('c', 1, size(C, 2)), '}', newline];
else
    s = [s, columnsSpec, '}', newline];
end

for i = 1:size(C, 1)
    for j = 1:size(C, 2)
        s = [s, C{i, j}];
        if j < size(C, 2)
            s = [s, ' & '];
        elseif i < size(C, 1)
            s = [s, ' \\', newline];
        else
            s = [s, newline];
        end
    end
end
s = [s, '\end{tabular}'];

%% enregistrement

fid = fopen(texPath, 'wt');
fprintf(fid, '%s', s);
fclose(fid);











