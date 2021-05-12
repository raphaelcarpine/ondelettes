%% data

saveFigs = false;
saveFormats = {'fig', 'epsc', 'png'};
saveFormats = {'fig'};

displayLegend = true;

sizeFigs = 8/10 * [560, 420];
% sizeFigs = 5/7 * [560, 420];

% geometrie pont et position capteurs
L_p = 17.5; % longueur pont
l_p = 4.84; % largeur pont
pos_capt = [4.375, 4.84; 8.75, 4.84; 8.75, 0; 4.375, 0; 0, 0; -4.375, 0; 0, 4.84; -8.75, 4.84] + [8.75, 0];
realShapePlotPont = @(shape, figTitle) shapePlotPlate([L_p, l_p], pos_capt,...
    shape, figTitle, false); % deformees pont

%% data

% files
dataFolder = 'pont sens/article analyse donnees reelles/save';

filesA = dir(fullfile(dataFolder, '*saveA*.mat')); % passage train
filesB = dir(fullfile(dataFolder, '*saveB*.mat')); % après

saveFolder = 'pont sens/article analyse donnees reelles/figs';

% sort
for kt = 1:length(filesA)
    filesAn(kt) = str2double(filesA(kt).name(4:end-9));
end
[~, I] = sort(filesAn);
filesA = filesA(I);

for kt = 1:length(filesB)
    filesBn(kt) = str2double(filesB(kt).name(4:end-9));
end
[~, I] = sort(filesBn);
filesB = filesB(I);

% variables
Nch = 8;

NtrainsA = length(filesA);
NtrainsB = length(filesB);

F0b = readtable(fullfile('pont sens/article analyse donnees reelles/FnQ_B.xlsx'));
NmodesA = 1;
NmodesB = 3;

ShapesA = nan(NtrainsA, NmodesA, Nch);
ShapesB = nan(NtrainsB, NmodesB, Nch);
GraphsA = cell(NtrainsA, 2);
GraphsB = cell(NtrainsB, 2);
GraphsB2 = cell(NtrainsB, 2);

% loading
for kt = 1:NtrainsA
    load(fullfile(dataFolder, filesA(kt).name));
    
    ShapesA(kt, :, :) = shapes;
    GraphsA{kt, 1} = ampl_graph;
    GraphsA{kt, 2} = freq_graph;
end
for kt = 1:NtrainsB
    load(fullfile(dataFolder, filesB(kt).name));
    
    ShapesB(kt, :, :) = shapes(1:NmodesB, :);
    GraphsB{kt, 1} = ampl_graph;
    GraphsB{kt, 2} = freq_graph;
    GraphsB2{kt, 1} = time_graph2;
    GraphsB2{kt, 2} = ampl_graph2;
end

%% figures

amplScale = 'log';
fig = figure;
c = colororder;
delete(fig);
c = [c; 0.8819 0.1904 0.4607; 0.6692 0.3689 0.9816];

for k = 1:min(NmodesA, NmodesB)
    fig0(k) = figure('Name', ['mode ', num2str(k), ' (pendant et après)']);
    fig0(k).Position(3:4) = sizeFigs;
    colororder(fig0(k), c);
    ax0(k) = axes;
    hold(ax0(k), 'on');
    xlabel(ax0(k), 'Amplitude [m/s²]');
    ylabel(ax0(k), 'Frequency [Hz]');
    set(ax0(k), 'XScale', amplScale);
end
for k = 1:NmodesA
    figA(k) = figure('Name', ['mode ', num2str(k), ' (pendant)']);
    figA(k).Position(3:4) = sizeFigs;
    colororder(figA(k), c);
    axA(k) = axes;
    hold(axA(k), 'on');
    xlabel(axA(k), 'Amplitude [m/s²]');
    ylabel(axA(k), 'Frequency [Hz]');
    set(axA(k), 'XScale', amplScale);
end
for k = 1:NmodesB
    figB(k) = figure('Name', ['mode ', num2str(k), ' (après)']);
    figB(k).Position(3:4) = sizeFigs;
    colororder(figB(k), c);
    axB(k) = axes;
    hold(axB(k), 'on');
    xlabel(axB(k), 'Amplitude [m/s²]');
    ylabel(axB(k), 'Frequency [Hz]');
    set(axB(k), 'XScale', amplScale);
end
for k = 1:NmodesB
    figB2(k) = figure('Name', ['mode ', num2str(k), ' (après)']);
    figB2(k).Position(3:4) = sizeFigs;
    colororder(figB2(k), c);
    axB2(k) = axes;
    hold(axB2(k), 'on');
    xlabel(axB2(k), 'Time [s]');
    ylabel(axB2(k), 'Amplitude [m/s²]');
    set(axB2(k), 'YScale', amplScale);
end

%% plot ampl freq

% passage
for kt = 1:NtrainsA
    for km = 1:NmodesA
        plt0 = plot(axA(km), GraphsA{kt, 1}{km}, GraphsA{kt, 2}{km}, 'DisplayName', ['train ', num2str(kt)]);
        if any(isnan(GraphsA{kt, 1}{km}))
            plt0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        
        if km <= NmodesB
            plot(ax0(km), GraphsA{kt, 1}{km}, GraphsA{kt, 2}{km}, 'DisplayName', ['train ', num2str(kt)]);
        end
    end
end
% reset couleurs
for km = 1:min(NmodesA, NmodesB)
    ax0(km).ColorOrderIndex = 1;
end
% apres passage
for kt = 1:NtrainsB
    for km = 1:NmodesB
        plt0 = plot(axB(km), GraphsB{kt, 1}{km}, GraphsB{kt, 2}{km}, 'DisplayName', ['train ', num2str(kt)]);
        if any(isnan(GraphsB{kt, 1}{km}))
            plt0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        
        plt0 = plot(axB2(km), GraphsB2{kt, 1}{km} - GraphsB2{kt, 1}{km}(1),...
            GraphsB{kt, 1}{km}, 'DisplayName', ['train ', num2str(kt)]);
        if any(isnan(GraphsB2{kt, 1}{km}))
            plt0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        
        if km <= NmodesA
            plt0 = plot(ax0(km), GraphsB{kt, 1}{km}, GraphsB{kt, 2}{km});
            plt0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
end

if displayLegend
%     XLimLegend0 = [0.018];
%     XLimLegendA = [nan];
%     XLimLegendB = [nan, 0.027, 0.02];
    for km = 1:min(NmodesA, NmodesB)
        legend(ax0(km), 'FontSize', 8);
    end
    for km = 1:NmodesA
        legend(axA(km), 'FontSize', 8);
    end
    for km = 1:NmodesB
        legend(axB(km), 'FontSize', 8);
        legend(axB2(km), 'FontSize', 8);
    end
end



%% plot shape

% passage
for km = 1:NmodesA
    shape = zeros(1, 1, Nch);
    nshape = 0;
    for kt = 1:size(ShapesA, 1)
        if ~isnan(ShapesA(kt, km, 1))
            shape2 = ShapesA(kt, km, :);
            if sum(real(shape) .* real(shape2), 'all') < 0
                shape2 = -shape2;
            end
            shape = shape + shape2;
            nshape = nshape +1;
        end
    end
    shape = shape / nshape;
    if -min(real(shape)) > max(real(shape))
        shape = -shape;
    end
    figShapeA(km) = realShapePlotPont(real(shape), ['shapeA_mode', num2str(km)]);
    figShapeA(km).Position(3:4) = sizeFigs;
end
% après passage
for km = 1:NmodesB
    shape = zeros(1, 1, Nch);
    nshape = 0;
    for kt = 1:size(ShapesB, 1)
        if ~isnan(ShapesB(kt, km, 1))
            shape2 = ShapesB(kt, km, :);
            if sum(real(shape) .* real(shape2), 'all') < 0
                shape2 = -shape2;
            end
            shape = shape + shape2;
            nshape = nshape +1;
        end
    end
    shape = shape / nshape;
    if -min(real(shape)) > max(real(shape))
        shape = -shape;
    end
    figShapeB(km) = realShapePlotPont(real(shape), ['shapeB_mode', num2str(km)]);
    figShapeB(km).Position(3:4) = sizeFigs;
end


%% save

if saveFigs
    for kf = 1:length(saveFormats)
        % graphs ampl freq
        for km = 1:NmodesA
            saveas(figA(km), fullfile(saveFolder, ['AmplFreqA_mode', num2str(km)]), saveFormats{kf});
        end
        for km = 1:NmodesB
            saveas(figB(km), fullfile(saveFolder, ['AmplFreqB_mode', num2str(km)]), saveFormats{kf});
            saveas(figB2(km), fullfile(saveFolder, ['TimeAmplB_mode', num2str(km)]), saveFormats{kf});
        end
        for km = 1:min(NmodesA, NmodesB)
            saveas(fig0(km), fullfile(saveFolder, ['AmplFreqAB_mode', num2str(km)]), saveFormats{kf});
        end
        % shapes
        for km = 1:NmodesA
            saveas(figShapeA(km), fullfile(saveFolder, ['shapeA_mode', num2str(km)]), saveFormats{kf});
        end
        for km = 1:NmodesB
            saveas(figShapeB(km), fullfile(saveFolder, ['shapeB_mode', num2str(km)]), saveFormats{kf});
        end
    end
end
















