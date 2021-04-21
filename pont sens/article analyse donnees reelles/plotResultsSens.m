%% data

saveFigs = false;
saveFormats = {'fig', 'eps', 'png'};

% geometrie pont et position capteurs
L_p = 17.5; % longueur pont
l_p = 4.84; % largeur pont
pos_capt = [4.375, 0; 8.75, 0; 8.75, 4.84; 4.375, 4.84; 0, 4.84; -4.375, 4.84; 0, 0; -8.75, 0] + [8.75, 0];
realShapePlotPont = @(shape, figTitle) shapePlotPlate([L_p, l_p], pos_capt,...
    shape * sign(real(shape(5))), figTitle, false); % deformees pont

%% data

% files
dataFolder = 'pont sens/article analyse donnees reelles/save';

filesA = dir(fullfile(dataFolder, '*saveA*.mat')); % passage train
filesB = dir(fullfile(dataFolder, '*saveB*.mat')); % après

saveFolder = 'pont sens/article analyse donnees reelles/figs';

% variables
Nch = 8;

NtrainsA = length(filesA);
NtrainsB = length(filesB);

F0b = readtable(fullfile('pont sens/article analyse donnees reelles/FnQ_B.xlsx'));
NmodesA = 1;
NmodesB = length(F0b.f0);

ShapesA = nan(NtrainsA, NmodesA, Nch);
ShapesB = nan(NtrainsB, NmodesB, Nch);
GraphsA = cell(NtrainsA, 2);
GraphsB = cell(NtrainsB, 2);

% loading
for kt = 1:NtrainsA
    load(fullfile(dataFolder, filesA(kt).name));
    
    ShapesA(kt, :, :) = shapes;
    GraphsA{kt, 1} = ampl_graph;
    GraphsA{kt, 2} = freq_graph;
end
for kt = 1:NtrainsB
    load(fullfile(dataFolder, filesB(kt).name));
    
    ShapesB(kt, :, :) = shapes;
    GraphsB{kt, 1} = ampl_graph;
    GraphsB{kt, 2} = freq_graph;
end

%% figures

amplScale = 'log';

for k = 1:min(NmodesA, NmodesB)
    fig0(k) = figure('Name', ['mode ', num2str(k), ' (pendant et après)']);
    ax0(k) = axes;
    hold(ax0(k), 'on');
    xlabel(ax0(k), 'Amplitude [m/s²]');
    ylabel(ax0(k), 'Frequency [Hz]');
    set(ax0(k), 'XScale', amplScale);
end
for k = 1:NmodesA
    figA(k) = figure('Name', ['mode ', num2str(k), ' (pendant)']);
    axA(k) = axes;
    hold(axA(k), 'on');
    xlabel(axA(k), 'Amplitude [m/s²]');
    ylabel(axA(k), 'Frequency [Hz]');
    set(axA(k), 'XScale', amplScale);
end
for k = 1:NmodesB
    figB(k) = figure('Name', ['mode ', num2str(k), ' (après)']);
    axB(k) = axes;
    hold(axB(k), 'on');
    xlabel(axB(k), 'Amplitude [m/s²]');
    ylabel(axB(k), 'Frequency [Hz]');
    set(axB(k), 'XScale', amplScale);
end

%% plot ampl freq

% passage
for kt = 1:NtrainsA
    for km = 1:NmodesA
        plot(axA(km), GraphsA{kt, 1}{km}, GraphsA{kt, 2}{km});
        if km <= NmodesB
            plot(ax0(km), GraphsA{kt, 1}{km}, GraphsA{kt, 2}{km});
        end
    end
end
% apres passage
for kt = 1:NtrainsB
    for km = 1:NmodesB
        plot(axB(km), GraphsB{kt, 1}{km}, GraphsB{kt, 2}{km});
        if km <= NmodesA
            plot(ax0(km), GraphsB{kt, 1}{km}, GraphsB{kt, 2}{km});
        end
    end
end



%% plot shape

% passage
for km = 1:NmodesA
    shape = zeros(1, 1, Nch);
    nshape = 0;
    for kt = 1:size(ShapesA, 1)
        if ~isnan(ShapesA(kt, km, 1))
            shape = shape + ShapesA(kt, km, :) * sign(real(ShapesA(kt, km, 5)));
            nshape = nshape +1;
        end
    end
    shape = shape / nshape;
    figShapeA(km) = realShapePlotPont(real(shape), ['shapeA_mode', num2str(km)]);
end
% après passage
for km = 1:NmodesB
    shape = zeros(1, 1, Nch);
    nshape = 0;
    for kt = 1:size(ShapesB, 1)
        if ~isnan(ShapesB(kt, km, 1))
            shape = shape + ShapesB(kt, km, :) * sign(real(ShapesB(kt, km, 5)));
            nshape = nshape +1;
        end
    end
    shape = shape / nshape;
    figShapeB(km) = realShapePlotPont(real(shape), ['shapeB_mode', num2str(km)]);
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
        end
        for km = 1:min(NmodesA, NmodesB)
            saveas(fig0(km), fullfile(saveFolder, ['AmplFreqAB_mode', num2str(km)]), saveFormats{kf});
        end
        % shapes
        for km = 1:length(figShapeA)
            saveas(figShapeA(km), fullfile(saveFolder, ['shapeA_mode', num2str(km)]), saveFormats{kf});
        end
        for km = 1:length(figShapeB)
            saveas(figShapeB(km), fullfile(saveFolder, ['shapeB_mode', num2str(km)]), saveFormats{kf});
        end
    end
end
















