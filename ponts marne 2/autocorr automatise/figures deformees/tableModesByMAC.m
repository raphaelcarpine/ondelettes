%% chargement data

dataPath = 'ponts marne 2\autocorr automatise\ondelette\save';

filesData = dir(fullfile(dataPath, '*.mat'));
filesPaths = {filesData.folder};
filesNames = {filesData.name};
for kfile = 1:length(filesNames)
    filesPaths{kfile} = fullfile(filesPaths{kfile}, filesNames{kfile});
    [~, filesNames{kfile}] = fileparts(filesNames{kfile});
end

ShapesCell = cell(1, length(filesNames));
FreqsCell = cell(1, length(filesNames));
DampsCell = cell(1, length(filesNames));

for kfile = 1:length(filesNames)
    load(filesPaths{kfile});
    ShapesCell{kfile} = Shapes;
    FreqsCell{kfile} = Freqs;
    DampsCell{kfile} = Damps;
end


%% calcul tableaux

modesMat = groupShapesByMac(ShapesCell);

modesFreqsMat = nan(size(modesMat));
modesDampsMat = nan(size(modesMat));
for kfile = 1:size(modesMat, 2)
    for kfreq = 1:size(modesMat, 1)
        if ~isnan(modesMat(kfreq, kfile))
            modesFreqsMat(kfreq, kfile) = FreqsCell{kfile}(modesMat(kfreq, kfile));
            modesDampsMat(kfreq, kfile) = DampsCell{kfile}(modesMat(kfreq, kfile));
        end
    end
end

% tri par freq croissante
freqsMoy = mean(modesFreqsMat, 2, 'omitnan');
[~, If] = sort(freqsMoy);
modesMat = modesMat(If, :);
modesFreqsMat = modesFreqsMat(If, :);
modesDampsMat = modesDampsMat(If, :);

% tableaux
modesTable = array2table(modesMat);
modesTable.Properties.VariableNames = filesNames;
modesFreqsTable = array2table(modesFreqsMat);
modesFreqsTable.Properties.VariableNames = filesNames;
modesDampsTable = array2table(modesDampsMat);
modesDampsTable.Properties.VariableNames = filesNames;


%% preparation tableau latex deformees

saveFolderLatex = 'ponts marne 2\autocorr automatise\figures deformees';
saveFolderFigsName = 'figures';
saveFolderFigs = fullfile(saveFolderLatex, saveFolderFigsName);

dimensionsShapes2; % script deformees modales

% enregistrement figure
for kmod = 1:size(modesMat, 1)
    for kfile = 1:size(modesMat, 2)
        if isnan(modesMat(kmod, kfile))
            continue
        end
        freq = modesFreqsMat(kmod, kfile);
        damp = modesDampsMat(kmod, kfile);
        shape = ShapesCell{kfile}(:, modesMat(kmod, kfile));
        In = nonPropIndex(shape);
        
%         shape = shape * sign(real(shape(8)));
        
        figName = [filesNames{kfile}, sprintf('_mode%02d', kmod)];
        fig = shapePlotBridge(real(shape), figName);
        set(fig, 'CloseRequestFcn', @closeFct);
        close(fig);
    end
end

%%

% creation fichier tex
xWidth = 1/(size(modesFreqsMat, 1)+1);

s = '\begin{tabular}{';
s = [s, 'l', repmat('|c', 1, size(modesFreqsMat, 1)), '}', newline];

for kmod = 1:size(modesFreqsMat, 1)
    s = [s, ' & ', num2str(kmod)];
end
s = [s, ' \\ \hline', newline];

for kfile = 1:size(modesFreqsMat, 2)
    titreLigne = filesNames{kfile};
    titreLigne = replace(titreLigne, 'modes_', '');
    titreLigne = replace(titreLigne, '_', ' ');
    s = [s, '\multirow{2}{*}{\rotatebox[origin=c]{90}{', titreLigne, '}} & '];
    for kmod = 1:size(modesFreqsMat, 1)
        if ~isnan(modesFreqsMat(kmod, kfile))
            freq = modesFreqsMat(kmod, kfile);
            damp = modesDampsMat(kmod, kfile);
            shape = ShapesCell{kfile}(:, modesMat(kmod, kfile));
            In = nonPropIndex(shape);
            
            s = [s, sprintf('f=%.2fHz, z=%.2f%%, I=%.2f%%', [freq, 100*damp, 100*In])];
        end
        
        if kmod < size(modesFreqsMat, 1)
            s = [s, ' & '];
        else
            s = [s, ' \\', newline];
        end
    end
    
    s = [s, ' & '];
    for kmod = 1:size(modesFreqsMat, 1)
        if ~isnan(modesFreqsMat(kmod, kfile))
            figFilePath = [filesNames{kfile}, sprintf('_mode%02d.png', kmod)];
            figFilePath = [saveFolderFigsName, '/', figFilePath];
            s = [s, '\includegraphics[width=', num2str(xWidth), '\linewidth]{', figFilePath, '}'];
        end
            
        
        if kmod < size(modesFreqsMat, 1)
            s = [s, ' & '];
        else
            s = [s, ' \\ \hline', newline];
        end
    end
end
s = [s, '\end{tabular}'];
s = replace(s, '%', '\%');

% ecriture ds le fichier
fileID = fopen(fullfile(saveFolderLatex, 'docModes.tex'), 'r');
s0 = fscanf(fileID, '%c');
ki = strfind(s0, '\begin{tabular}');
kf = strfind(s0, '\end{tabular}');
fclose(fileID);
try
    if length(ki) ~= 1 || length(kf) ~= 1
        error(' ');
    end

    S = [s0(1:ki-1), s, s0(kf+13:end)];
    fileID = fopen(fullfile(saveFolderLatex, 'docModes.tex'), 'w');
    fprintf(fileID, '%s', S);
    fclose(fileID);
catch
    fprintf(2, 'problème lecture/écriture fichier\n');
end






%%

function closeFct(fig, ~)
    saveFolderLatex = 'ponts marne 2\autocorr automatise\figures deformees';
    saveFolderFigsName = 'figures';
    saveFolderFigs = fullfile(saveFolderLatex, saveFolderFigsName);
    saveas(fig, fullfile(saveFolderFigs, fig.Name), 'png');
    delete(fig);
end

