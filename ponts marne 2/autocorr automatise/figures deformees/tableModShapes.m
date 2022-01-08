saveFigs = 0; % réenregistrement des figures


%% chargement data

[filesNames, FreqsCell, DampsCell, ShapesCell] = getCWTResults();


%% calcul température moyenne

TempMoyCell = cell(size(filesNames));

for bridge = 1:10
    % chargement data pont pour heures début - fin
    dataFilePath = choixData(bridge);
    load(dataFilePath);
    disp(startDate);
    X = X.';
    T = T.';
    [X, T] = removeRedundantData(X, T);
    [~, T] = removeNanSignal(X, T);
    [~, dataFileName] = fileparts(dataFilePath);
    disp(' ');
    
    % chargement température
    [timeTemp, Temp] = getTemperature(startDate);
%     figure;
%     plot(timeTemp, Temp);
    
    % calcul moyenne température
    if timeTemp(1) > startDate + T/(24*3600)
        warning('timeTemp(1) > T(1)');
    end
    Tmoy = interp1(timeTemp, Temp, startDate + T/(24*3600));
%     figure;
%     plot(startDate + T/(24*3600), Tmoy);
%     xlim(startDate + T([1, end])/(24*3600));
    Tmoy = mean(Tmoy, 'omitnan');
    
    % recherche fichier correspondant
    kfile = 1;
    while kfile <= length(filesNames) && ~strcmp(filesNames{kfile}(7:end), dataFileName)
        kfile = kfile + 1;
    end
    if kfile > length(filesNames)
        warning('fichier non trouvé');
    else
        TempMoyCell{kfile} = Tmoy;
    end
end


%% calcul tableaux

modesMat = groupModShapes(FreqsCell, ShapesCell);

modesFreqsMat = nan(size(modesMat));
modesDampsMat = nan(size(modesMat));
for kfreq = 1:size(modesMat, 1)
    for kfile = 1:size(modesMat, 2)
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


%% preparation tableau latex et figures deformees

saveFolderLatex = 'ponts marne 2\autocorr automatise\figures deformees';
saveFolderFigsName = 'figures';
saveFolderFigs = fullfile(saveFolderLatex, saveFolderFigsName);

dimensionsShapes2; % script deformees modales

% enregistrement figure
if saveFigs
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


%% graphs freqs & amorts

grpPonts = {[1], [2 3], [4 5 6], [7 8], [9]};
MarkersOrdre = {'+', 'x', '*'};

filesNames2 = cellfun(@(s) replace(replace(s, 'modes_', ''), '_', ' '), filesNames, 'UniformOutput', false);
filesNames2Temp = filesNames2;
for k = 1:length(filesNames2Temp)
    filesNames2Temp{k} = [filesNames2Temp{k}, sprintf(', %.1f°C', TempMoyCell{k})];
end

modesFreqsDevMat = modesFreqsMat./mean(modesFreqsMat, 2, 'omitnan') - 1;
modesDampsDevMat = modesDampsMat./mean(modesDampsMat, 2, 'omitnan') - 1;

figure;
for kpont = 1:length(grpPonts)
    for kf = 1:length(grpPonts{kpont})
        kfile = grpPonts{kpont}(kf);
        set(gca, 'ColorOrderIndex', kpont)
        plot(100*modesFreqsDevMat(:, kfile), 'Marker', MarkersOrdre{kf}, 'MarkerSize', 8);
        hold on
    end
end
xlabel('mode');
ylabel('Df [%]');
xlim([2.5 10.5]);
legend(filesNames2Temp);
selectLine();

figure;
for kpont = 1:length(grpPonts)
    for kf = 1:length(grpPonts{kpont})
        kfile = grpPonts{kpont}(kf);
        set(gca, 'ColorOrderIndex', kpont)
        plot(100*modesDampsDevMat(:, kfile), 'Marker', MarkersOrdre{kf}, 'MarkerSize', 8);
        hold on
    end
end
xlabel('mode');
ylabel('Dz [%]');
xlim([2.5 10.5]);
legend(filesNames2Temp);
selectLine();

%% dépendance frequence et amort en température

dfdT_mod = nan(1, 10);
f0_mod = nan(1, 10);
dzdT_mod = nan(1, 10);
z0_mod = nan(1, 10);

% reg lin par mode
plotRegLin = false;
for kmode = 3:10
    if plotRegLin
        figure('Name', sprintf('Dépendance fréquence propore en température, mode %d', kmode));
        axF = axes();
        hold(axF, 'on');
        figure('Name', sprintf('Dépendance amortissement en température, mode %d', kmode));
        axZ = axes();
        hold(axZ, 'on');
        for kpont = 1:length(grpPonts)
            %             if length(grpPonts{kpont}) <= 1
            %                 continue
            %             end
            %         Tmoy = mean([TempMoyCell{grpPonts{kpont}}]);
            for kfile = grpPonts{kpont}
                set(axF, 'ColorOrderIndex', kpont)
                plot(axF, TempMoyCell{kfile}*ones(size(modesFreqsDevMat(kmode, kfile))), 100*modesFreqsDevMat(kmode, kfile), '+');
                set(axZ, 'ColorOrderIndex', kpont)
                plot(axZ, TempMoyCell{kfile}*ones(size(modesDampsDevMat(kmode, kfile))), 100*modesDampsDevMat(kmode, kfile), '+');
            end
        end
    end
    
    % reg lin
    DfreqMode = [];
    DdampMode = [];
    tempMode = [];
    for kpont = 1:length(grpPonts)
        for kfile = grpPonts{kpont}
            if ~isnan(modesFreqsDevMat(kmode, kfile))
                tempMode(end+1) = TempMoyCell{kfile};
                DfreqMode(end+1) = modesFreqsDevMat(kmode, kfile);
                DdampMode(end+1) = modesDampsDevMat(kmode, kfile);
            end
        end
    end
    coeffs = [ones(size(tempMode)); tempMode].' \ DfreqMode.';
    f0 = coeffs(1);
    dfdT = coeffs(2);
    coeffs = [ones(size(tempMode)); tempMode].' \ DdampMode.';
    z0 = coeffs(1);
    dzdT = coeffs(2);
    if plotRegLin
        TempLim = get(axF, 'xlim');
        plot(axF, TempLim, 100*(f0 + dfdT*TempLim), '--r');
        TempLim = get(axZ, 'xlim');
        plot(TempLim, 100*(z0 + dzdT*TempLim), '--r');
    end
    
    dfdT_mod(kmode) = dfdT;
    f0_mod(kmode) = f0;
    dzdT_mod(kmode) = dzdT;
    z0_mod(kmode) = z0;
end

% pente moyenne
dfdT_moy = mean(dfdT_mod, 'omitnan');
dzdT_moy = mean(dzdT_mod, 'omitnan');

figure;
plot(100*dfdT_mod, '+');
% yline(100*dfdT_moy, '--r');
xlabel('mode');
ylabel('(df_n/dT)/f_n [%/°C]');
xlim([2.5, 10.5]);
ax = gca;
ax.YLim(2) = 0;

figure;
plot(100*dzdT_mod, '+');
% yline(100*dzdT_moy, '--r');
xlabel('mode');
ylabel('(d\zeta_n/dT)/\zeta_n [%/°C]');
xlim([2.5, 10.5]);

%% recalage température

T0_mod = - f0_mod ./ dfdT_mod;
modesFreqsDevMatCorr = modesFreqsDevMat;
for kmode = 1:10
    T0 = T0_mod(kmode);
    dfdT = dfdT_moy;
%     dfdT = dfdT_mod(kmode);
    modesFreqsDevMatCorr(kmode, :) = modesFreqsDevMatCorr(kmode, :) - dfdT * ([TempMoyCell{:}] - T0);
end

figure;
for kpont = 1:length(grpPonts)
    for kf = 1:length(grpPonts{kpont})
        kfile = grpPonts{kpont}(kf);
        set(gca, 'ColorOrderIndex', kpont)
        plot(100*modesFreqsDevMatCorr(:, kfile), 'Marker', MarkersOrdre{kf}, 'MarkerSize', 8);
        hold on
    end
end
xlabel('mode');
ylabel('Df [%]');
xlim([2.5 10.5]);
legend(filesNames2Temp);
selectLine();


%% moyennes

dimensionsShapes2

for kmode = 1:10
    f_moy = 0;
    z_moy = 0;
    shape_moy = zeros(size(ShapesCell{1}, 1), 1);
    n_moy = 0;
    for kfile = 1:length(FreqsCell)
        kmode2 = modesMat(kmode, kfile);
        if ~isnan(kmode2) && ~isnan(FreqsCell{kfile}(kmode2))
            f_moy = f_moy + FreqsCell{kfile}(kmode2);
            z_moy = z_moy + DampsCell{kfile}(kmode2);
            shape_moy = shape_moy + ShapesCell{kfile}(:, kmode2);
            n_moy = n_moy + 1;
        end
    end
    f_moy = f_moy / n_moy;
    z_moy = z_moy / n_moy;
    shape_moy = shape_moy / n_moy;
    In = nonPropIndex(shape_moy);
    figName = sprintf('deformée moyenne mode %d, f=%.2fHz, z=%.1f%%, I=%.1f%%', [kmode, f_moy, 100*z_moy, 100*In]);
    shapePlotBridge(real(shape_moy), figName);
end


%%

function closeFct(fig, ~)
    saveFolderLatex = 'ponts marne 2\autocorr automatise\figures deformees';
    saveFolderFigsName = 'figures';
    saveFolderFigs = fullfile(saveFolderLatex, saveFolderFigsName);
    saveas(fig, fullfile(saveFolderFigs, fig.Name), 'png');
    delete(fig);
end

