%% data

saveFigs = false;
saveFormats = {'fig', 'epsc', 'png'};
saveFormats = {'fig'};

displayLegend = true;

sizeFigs = 8/10 * [560, 420];
% sizeFigs = 5/7 * [560, 420];

dt = 1/4096;

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
GraphsA1 = cell(NtrainsA, 2);
GraphsB1 = cell(NtrainsB, 2);
GraphsB2 = cell(NtrainsB, 2);

% loading
for kt = 1:NtrainsA
    load(fullfile(dataFolder, filesA(kt).name));
    
    ShapesA(kt, :, :) = shapes;
    GraphsA{kt, 1} = ampl_graph;
    GraphsA{kt, 2} = freq_graph;
    
    % calcul amort
    GraphsA1{kt, 1} = cell(1, length(ampl_graph));
    GraphsA1{kt, 2} = cell(1, length(ampl_graph));
    for km = 1:length(ampl_graph)
        GraphsA1{kt, 1}{km} = ampl_graph{km};
        omega = 2*pi*freq_graph{km};
        if length(omega) >= 2
            logA = log(ampl_graph{km});
            dlogA = diff(logA);
            dlogA = [dlogA(1), 1/2*(dlogA(1:end-1)+dlogA(2:end)), dlogA(end)];
            dlogA = dlogA / dt;
            z = - dlogA ./ abs(dlogA + 1i*omega);
            % lissage
            deltaT = 0.1;
            Nhalf = ceil(5*deltaT/dt);
            Nhalf = min(Nhalf, length(z));
            filterCoeffs = nan(1, 2*Nhalf+1);
            for kfilter = 1:length(filterCoeffs)
                filterCoeffs(kfilter) = exp( -((kfilter-1-Nhalf)*dt)^2 / (2*deltaT^2));
            end
            filterCoeffs = filterCoeffs / sum(filterCoeffs);
            dampingData = z;
            dampingData = [dampingData(1)*ones(1, Nhalf), dampingData, dampingData(end)*ones(1, Nhalf)];
            dampingData = filter(filterCoeffs, 1, dampingData);
            z = dampingData(2*Nhalf+1:end);
        elseif length(omega) == 1
            z= nan;
        else
            z = [];
        end
        GraphsA1{kt, 2}{km} = z;
    end
end
for kt = 1:NtrainsB
    load(fullfile(dataFolder, filesB(kt).name));
    
    ShapesB(kt, :, :) = shapes(1:NmodesB, :);
    GraphsB{kt, 1} = ampl_graph;
    GraphsB{kt, 2} = freq_graph;
    GraphsB2{kt, 1} = time_graph2;
    GraphsB2{kt, 2} = ampl_graph2;
    
    % calcul amort
    GraphsB1{kt, 1} = cell(1, length(ampl_graph));
    GraphsB1{kt, 2} = cell(1, length(ampl_graph));
    for km = 1:length(ampl_graph)
        GraphsB1{kt, 1}{km} = ampl_graph{km};
        omega = 2*pi*freq_graph{km};
        if length(omega) >= 2
            logA = log(ampl_graph{km});
            dlogA = diff(logA);
            dlogA = [dlogA(1), 1/2*(dlogA(1:end-1)+dlogA(2:end)), dlogA(end)];
            dlogA = dlogA / dt;
            z = - dlogA ./ abs(dlogA + 1i*omega);
            % lissage
            deltaT = 0.1;
            Nhalf = ceil(5*deltaT/dt);
            Nhalf = min(Nhalf, length(z));
            filterCoeffs = nan(1, 2*Nhalf+1);
            for kfilter = 1:length(filterCoeffs)
                filterCoeffs(kfilter) = exp( -((kfilter-1-Nhalf)*dt)^2 / (2*deltaT^2));
            end
            filterCoeffs = filterCoeffs / sum(filterCoeffs);
            dampingData = z;
            dampingData = [dampingData(1)*ones(1, Nhalf), dampingData, dampingData(end)*ones(1, Nhalf)];
            dampingData = filter(filterCoeffs, 1, dampingData);
            z = dampingData(2*Nhalf+1:end);
        elseif length(omega) == 1
            z = nan;
        else
            z = [];
        end
        GraphsB1{kt, 2}{km} = z;
    end
end

%% figures

amplScale = 'log';
fig = figure;
c = colororder;
delete(fig);
c = [c; 0.8819 0.1904 0.4607; 0.6692 0.3689 0.9816];

for k = 1:min(NmodesA, NmodesB) % ampl freq avant et après
    fig0(k) = figure('Name', ['mode ', num2str(k), ' (pendant et après)']);
    fig0(k).Position(3:4) = sizeFigs;
    colororder(fig0(k), c);
    ax0(k) = axes('Box', true);
    hold(ax0(k), 'on');
    xlabel(ax0(k), 'Amplitude [m/s²]');
    ylabel(ax0(k), 'Frequency [Hz]');
    set(ax0(k), 'XScale', amplScale);
    set(fig0(k), 'renderer', 'painters');
end
for k = 1:NmodesA % ampl freq avant
    figA(k) = figure('Name', ['mode ', num2str(k), ' (pendant)']);
    figA(k).Position(3:4) = sizeFigs;
    colororder(figA(k), c);
    axA(k) = axes('Box', true);
    hold(axA(k), 'on');
    xlabel(axA(k), 'Amplitude [m/s²]');
    ylabel(axA(k), 'Frequency [Hz]');
    set(axA(k), 'XScale', amplScale);
    set(figA(k), 'renderer', 'painters');
end
for k = 1:NmodesB % ampl freq après
    figB(k) = figure('Name', ['mode ', num2str(k), ' (après)']);
    figB(k).Position(3:4) = sizeFigs;
    colororder(figB(k), c);
    axB(k) = axes('Box', true);
    hold(axB(k), 'on');
    xlabel(axB(k), 'Amplitude [m/s²]');
    ylabel(axB(k), 'Frequency [Hz]');
    set(axB(k), 'XScale', amplScale);
    set(figB(k), 'renderer', 'painters');
end
for k = 1:min(NmodesA, NmodesB) % ampl damp avant et après
    fig01(k) = figure('Name', ['mode ', num2str(k), ' (pendant et après)']);
    fig01(k).Position(3:4) = sizeFigs;
    colororder(fig01(k), c);
    ax01(k) = axes('Box', true);
    hold(ax01(k), 'on');
    xlabel(ax01(k), 'Amplitude [m/s²]');
    ylabel(ax01(k), 'Damping ratio [%]');
    set(ax01(k), 'XScale', amplScale);
    set(fig01(k), 'renderer', 'painters');
end
for k = 1:NmodesA % ampl damp avant
    figA1(k) = figure('Name', ['mode ', num2str(k), ' (pendant)']);
    figA1(k).Position(3:4) = sizeFigs;
    colororder(figA1(k), c);
    axA1(k) = axes('Box', true);
    hold(axA1(k), 'on');
    xlabel(axA1(k), 'Amplitude [m/s²]');
    ylabel(axA1(k), 'Damping ratio [%]');
    set(axA1(k), 'XScale', amplScale);
    set(figA1(k), 'renderer', 'painters');
end
for k = 1:NmodesB % ampl damp après
    figB1(k) = figure('Name', ['mode ', num2str(k), ' (après)']);
    figB1(k).Position(3:4) = sizeFigs;
    colororder(figB1(k), c);
    axB1(k) = axes('Box', true);
    hold(axB1(k), 'on');
    xlabel(axB1(k), 'Amplitude [m/s²]');
    ylabel(axB1(k), 'Damping ratio [%]');
    set(axB1(k), 'XScale', amplScale);
    set(figB1(k), 'renderer', 'painters');
end
for k = 1:NmodesB % time ampl après
    figB2(k) = figure('Name', ['mode ', num2str(k), ' (après)']);
    figB2(k).Position(3:4) = sizeFigs;
    colororder(figB2(k), c);
    axB2(k) = axes('Box', true);
    hold(axB2(k), 'on');
    xlabel(axB2(k), 'Time [s]');
    ylabel(axB2(k), 'Amplitude [m/s²]');
    set(axB2(k), 'YScale', amplScale);
    set(figB2(k), 'renderer', 'painters');
end

%% plot ampl freq

% passage
for kt = 1:NtrainsA
    for km = 1:NmodesA
        % ampl freq
        plt0 = plot(axA(km), GraphsA{kt, 1}{km}, GraphsA{kt, 2}{km}, 'DisplayName', ['train ', num2str(kt)]);
        if any(isnan(GraphsA{kt, 1}{km}))
            plt0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        % ampl damp
        plt0 = plot(axA1(km), GraphsA1{kt, 1}{km}, 100*GraphsA1{kt, 2}{km}, 'DisplayName', ['train ', num2str(kt)]);
        if any(isnan(GraphsA1{kt, 1}{km}))
            plt0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        
        if km <= NmodesB
            % ampl freq
            plot(ax0(km), GraphsA{kt, 1}{km}, GraphsA{kt, 2}{km}, 'DisplayName', ['train ', num2str(kt)]);
            % ampl damp
            plot(ax01(km), GraphsA1{kt, 1}{km}, 100*GraphsA1{kt, 2}{km}, 'DisplayName', ['train ', num2str(kt)]);
        end
    end
end
% reset couleurs
for km = 1:min(NmodesA, NmodesB)
    ax0(km).ColorOrderIndex = 1;
    ax01(km).ColorOrderIndex = 1;
end
% apres passage
for kt = 1:NtrainsB
    for km = 1:NmodesB
        % ampl freq
        plt0 = plot(axB(km), GraphsB{kt, 1}{km}, GraphsB{kt, 2}{km}, 'DisplayName', ['train ', num2str(kt)]);
        if any(isnan(GraphsB{kt, 1}{km}))
            plt0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        % ampl freq
        plt0 = plot(axB1(km), GraphsB1{kt, 1}{km}, 100*GraphsB1{kt, 2}{km}, 'DisplayName', ['train ', num2str(kt)]);
        if any(isnan(GraphsB1{kt, 1}{km}))
            plt0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        % time ampl
        plt0 = plot(axB2(km), GraphsB2{kt, 1}{km} - GraphsB2{kt, 1}{km}(1),...
            GraphsB{kt, 1}{km}, 'DisplayName', ['train ', num2str(kt)]);
        if any(isnan(GraphsB2{kt, 1}{km}))
            plt0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        
        if km <= NmodesA
            plt0 = plot(ax0(km), GraphsB{kt, 1}{km}, GraphsB{kt, 2}{km}); % ampl freq
            plt0.Annotation.LegendInformation.IconDisplayStyle = 'off';
            plt0 = plot(ax01(km), GraphsB1{kt, 1}{km}, 100*GraphsB1{kt, 2}{km}); % ampl damp
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
        legend(ax01(km), 'FontSize', 8);
    end
    for km = 1:NmodesA
        legend(axA(km), 'FontSize', 8);
        legend(axA1(km), 'FontSize', 8);
    end
    for km = 1:NmodesB
        legend(axB(km), 'FontSize', 8);
        legend(axB1(km), 'FontSize', 8);
        legend(axB2(km), 'FontSize', 8);
    end
end

% limites
for k = 1:min(NmodesA, NmodesB) % ampl damp avant et après
    ax01(k).YLim(1) = 0;
end
for k = 1:NmodesA % ampl damp avant
    axA1(k).YLim(1) = 0;
end
for k = 1:NmodesB % ampl damp après
    axB1(k).YLim(1) = 0;
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
            set(figA(km), 'renderer', 'painters');
            saveas(figA(km), fullfile(saveFolder, ['AmplFreqA_mode', num2str(km)]), saveFormats{kf});
        end
        for km = 1:NmodesB
            set(figB(km), 'renderer', 'painters');
            saveas(figB(km), fullfile(saveFolder, ['AmplFreqB_mode', num2str(km)]), saveFormats{kf});
            set(figB2(km), 'renderer', 'painters');
            saveas(figB2(km), fullfile(saveFolder, ['TimeAmplB_mode', num2str(km)]), saveFormats{kf});
        end
        for km = 1:min(NmodesA, NmodesB)
            set(fig0(km), 'renderer', 'painters');
            saveas(fig0(km), fullfile(saveFolder, ['AmplFreqAB_mode', num2str(km)]), saveFormats{kf});
        end
        % graphs ampl damp
        for km = 1:NmodesA
            set(figA1(km), 'renderer', 'painters');
            saveas(figA1(km), fullfile(saveFolder, ['AmplDampA_mode', num2str(km)]), saveFormats{kf});
        end
        for km = 1:NmodesB
            set(figB1(km), 'renderer', 'painters');
            saveas(figB1(km), fullfile(saveFolder, ['AmplDampB_mode', num2str(km)]), saveFormats{kf});
        end
        for km = 1:min(NmodesA, NmodesB)
            set(fig01(km), 'renderer', 'painters');
            saveas(fig01(km), fullfile(saveFolder, ['AmplDampAB_mode', num2str(km)]), saveFormats{kf});
        end
        % shapes
        for km = 1:NmodesA
            set(figShapeA(km), 'renderer', 'painters');
            saveas(figShapeA(km), fullfile(saveFolder, ['shapeA_mode', num2str(km)]), saveFormats{kf});
        end
        for km = 1:NmodesB
            set(figShapeB(km), 'renderer', 'painters');
            saveas(figShapeB(km), fullfile(saveFolder, ['shapeB_mode', num2str(km)]), saveFormats{kf});
        end
    end
end
















