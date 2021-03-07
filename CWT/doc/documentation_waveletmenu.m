function documentation_waveletmenu(editMode)
%DOCUMENTATION Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
    editModeStr = questdlg('Mode', 'Normal mode or edit mode', 'Normal', 'Edit', 'Edit');
    switch editModeStr
        case 'Normal'
            editMode = false;
        case 'Edit'
            editMode = true;
    end
end

title = 'Documentation';
fig = figure('Name', title, 'MenuBar', 'none', 'numbertitle', 'off');

savedTag = true;

folderPath = 'CWT/doc/doc/';


%% arborescence

arbo = {
    {'Utilisation', 'uti', {
        {'Démarer', 'dem'},...
        {'Paramètres', 'prm'},...
        {'Extraction des données', 'dat'}
    }},...
    {'Documentation générale', 'gen', {
        {'Transformée en ondelettes', 'cwt'},...
        {'Effets de bord', 'edg'},...
        {'Calcul des ridges', 'rdg'},...
        {'Transformée de Fourier', 'fft'}
    }},...
    {'Interface', 'int', {
        {'parameters', 'prm'},...
        {'plots', 'plt'},...
        {'mode shapes', 'shp'}
    }},...
    {'Menus', 'men', {
        {'Options', 'opt', {
            {'Mother wavelet', 'wvl'},...
            {'Frequency scale', 'fsc'},...
            {'Set channels', 'sch'},...
            {'Set time limits', 'tli'},...
            {'Remove mean', 'rmv'},...
            {'Set ct', 'sct'},...
            {'Set cf', 'scf'}
        }},...
        {'Mode', 'mod', {
            {'Multi signal mode', 'mlt'},...
            {'Random decrement mode', 'rdm'},...
            {'Cross-corr mode', 'xcm'},...
            {'Cross-corr params', 'xcp', {
                {'Set max lag', 'smt'},...
                {'Bias', 'bia'},...
                {'Set nb of SV', 'nsv'},...
                {'SVD mode CWT', 'svw'},...
                {'SVD mode Fourier', 'svf'}
            }}
        }},...
        {'Ridges', 'rdg', {
            {'Frequency', 'frq'},...
            {'Phase', 'pha'},...
            {'Damping', 'dmp'},...
            {'Set time limits ridge', 'tli'},...
            {'Set ct ridge', 'ctr'},...
            {'Stop ridge when increasing', 'stp'},...
            {'Set ridge threshold', 'thr'},...
            {'Get average |CWT|', 'avg'},...
            {'Multiple axes', 'mlt'}
        }},...
        {'Fourier', 'fft', {
            {'Window', 'wnd'},...
            {'Averaging', 'avg'}
        }},...
        {'Tools', 'tol', {
            {'Bounds Q', ''},...
            {'Filtering', ''},...
            {'Regression', ''},...
            {'Plot extract', ''}
        }}
    }}
};


%% construction des tabs

    function createTabs(arbo, parent, parentTag) % création récursive des onglets dans les onglets
        tabgrp = uitabgroup(parent, 'Units', 'normalized', 'Position', [0 0 1 1]);
        for k = 1:length(arbo)
            tab = uitab(tabgrp, 'Title', arbo{k}{1}, 'UserData', arbo{k}{1}); % création de l'onglet
            tag = [parentTag, '_', arbo{k}{2}];
            if length(arbo{k}) == 2 % crétion de la page de doc
                if ~editMode % mode normal
                    try
                        str = fileread([folderPath, tag, '.txt']);
                    catch
                        str = 'Fichier source introuvable ou corrompu.';
                    end
                    annotation(tab, 'textbox','String', str, 'Units', 'normalized',...
                        'Position', [0 0 1 1], 'LineStyle', 'none', 'VerticalAlignment', 'top', 'margin', 10);
                else % mode edition
                    try
                        str = fileread(['CWT/doc/doc/', tag, '.txt']);
                    catch
                        str = '';
                    end
                    hiddenPan = annotation(tab, 'textbox','String', str, 'Units', 'normalized',...
                        'Position', [0 0 1 1], 'LineStyle', 'none', 'VerticalAlignment', 'top', 'margin', 10);
                    str_cell = {};
                    i_str = 1;
                    while i_str <= length(str)
                        if str(i_str) == newline
                            str_cell{end+1} = str(1:i_str-1);
                            str = str(i_str+1:end);
                            i_str = 1;
                        elseif double(str(i_str)) == 13 % cariage return
                            str = [str(1:i_str-1), str(i_str+1:end)];
                        else
                            i_str = i_str+1;
                        end
                    end
                    str_cell{end+1} = str;
                    uicontrol(tab, 'Style', 'edit', 'String', str_cell, 'Units', 'normalized', 'Position', [0 0 1 1],...
                        'FontSize', 10, 'Max', 2, 'HorizontalAlignment', 'left', 'Tag', tag, 'UserData', hiddenPan);
                end
            elseif length(arbo{k}) == 3 % création récursive des onglets dans les onglets
                createTabs(arbo{k}{3}, tab, tag)
            end
        end
    end

createTabs(arbo, fig, 'doc');


%% scroll

allTxtPans = findall(fig, 'Type', 'textboxshape');
allTabGrps = findall(fig, 'Type', 'uitabgroup');

    function tab = getCurrentTab(parentTab)
        if nargin == 0
            parentTab = fig;
        end
        grptab = findobj(parentTab, 'Type', 'uitabgroup', '-depth', 1);
        if isempty(grptab)
            tab = parentTab;
        elseif length(grptab) == 1
            tab = getCurrentTab(grptab.SelectedTab);
        else
            error('');
        end
    end


    function scrollCallback(~, callbackdata)
        currentTab = getCurrentTab();
        txtbox = findall(currentTab, 'Type', 'textboxshape');
        txtbox.Position(2) = max(txtbox.Position(2) + 0.1*callbackdata.VerticalScrollCount, 0);
        drawnow;
    end

    function scrollReset(~, ~)
        for i = 1:length(allTxtPans)
            allTxtPans(i).Position(2) = 0;
        end
        drawnow;
    end

fig.WindowScrollWheelFcn = @scrollCallback;
for k_tbgrp = 1:length(allTabGrps)
    allTabGrps(k_tbgrp).SelectionChangedFcn = @scrollReset;
end


%% save

allEditPans = findall(fig, 'Type', 'uicontrol', '-and', 'Style', 'edit');

    function toggleSaveMode(tag, editPan)
        if tag % saving
            for k = 1:length(allEditPans)
                Str = get(allEditPans(k), 'String');
                str = '';
                for k_line = 1:length(Str)
                    str = [str, Str{k_line}];
                    if k_line < size(Str, 1)
                        str = [str, newline];
                    end
                end
                editTag = get(allEditPans(k), 'Tag');
                fid = fopen([folderPath, editTag, '.txt'], 'wt');
                fprintf(fid, '%s', str);
                fclose(fid);
            end
        end
        
        savedTag = tag;
        
        if savedTag % name change (*)
            for k = 1:length(allEditPans)
                allEditPans(k).Parent.Title = allEditPans(k).Parent.UserData;
            end
            set(fig, 'Name', title);
        else
            editPan.Parent.Title = [editPan.Parent.UserData, '*'];
            set(fig, 'Name', [title, '*']);
        end
    end

for j = 1:length(allEditPans) % not saved
    allEditPans(j).KeyPressFcn = @(~, ~) toggleSaveMode(false, allEditPans(j));
end

% saving
cm = uicontextmenu(fig);
cmsave = uimenu(cm, 'Text', 'Save', 'Accelerator', 'S');
cmsave.MenuSelectedFcn = @(~, ~) toggleSaveMode(true);
for j = 1:length(allEditPans)
    allEditPans(j).ContextMenu = cm;
end


%% view render (edit mode)

renderingMode = false;

    function toggleRender()
        scrollReset();
        renderingMode = ~renderingMode;
        if renderingMode
            for k = 1:length(allEditPans)
                editPan = allEditPans(k);
                hidenPan = editPan.UserData;
                Str = get(editPan, 'String');
                str = '';
                for k_line = 1:length(Str)
                    str = [str, Str{k_line}];
                    if k_line < size(Str, 1)
                        str = [str, newline];
                    end
                end
                hidenPan.String = str;
                editPan.Visible = 'off';
            end
        else
            for k = 1:length(allEditPans)
                allEditPans(k).Visible = 'on';
            end
        end
    end

cmrender = uimenu(cm, 'Text', 'Rendering mode');
cmrender.MenuSelectedFcn = @(~, ~) toggleRender();

if editMode
    for j = 1:length(allEditPans)
        allEditPans(j).UserData.ContextMenu = cm;
    end
end


%% closing

    function closeFcn(varargin)
        if ~savedTag
            selection = questdlg('Save before closing?', 'Saving', 'Yes','No', 'Cancel', 'Yes');
            switch selection
                case 'Yes'
                    toggleSaveMode(true);
                case 'No'
                case 'Cancel'
                    return
            end
        end
        delete(fig);
    end

fig.CloseRequestFcn = @closeFcn;
        


end

