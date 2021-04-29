%% inputs : 2, 6, 8, 9, 13, 15
N = 1;
mode = 4; % 0: affichage, 1: t0, 1.5 : t0 affichage, 2: ft, 3: passage train, 4 :apres train
saveResults = false;

channels = [3 4 5 6 7 9 10 15];
% temps_avant = dataArray{:, 1};
% c_02_avant = dataArray{:, 2};
% a_p12_avant = dataArray{:, 3};
% a_p15_avant = dataArray{:, 4};
% a_p13_avant = dataArray{:, 5};
% a_p10_avant = dataArray{:, 6};
% a_p07_avant = dataArray{:, 7};
% c_08_avant = dataArray{:, 8};
% a_p04_avant = dataArray{:, 9};
% a_p09_avant = dataArray{:, 10};
% c_11_avant = dataArray{:, 11};
% c_12_avant = dataArray{:, 12};
% c_13_avant = dataArray{:, 13};
% c_14_avant = dataArray{:, 14};
% a_p03_avant = dataArray{:, 15};
% c_16_avant = dataArray{:, 16};
% c_17_avant = dataArray{:, 17};

%% data
% longueurs
L_w = 18.7; % longueur wagon
L_t = 200.19 - 2*(5.02-1.5); % longueur train
L_loco = 1.5 + 14 + 1.5; % longueur locomotive
L_p = 17.5; % longueur pont

% geometrie pont et position capteurs
l_p = 4.84; % largeur pont
pos_capt = [4.375, 0; 8.75, 0; 8.75, 4.84; 4.375, 4.84; 0, 4.84; -4.375, 4.84; 0, 0; -8.75, 0] + [8.75, 0];
realShapePlotPont = @(shape, figTitle) shapePlotPlate([L_p, l_p], pos_capt,...
    shape * sign(real(shape(5))), figTitle, false); % deformees pont

% fréquences
ft0 = 4.3; % estimation ft

%% paths
dataFolder = 'pont sens/donnees reelles/data';
dataName = ['TGV', num2str(N), 'A.mat'];
dataPath = fullfile(dataFolder, dataName);

T0Folder = 'pont sens/article analyse donnees reelles';
T0Name = 'T0.xlsx';
T0Path = fullfile(T0Folder, T0Name);

FtFolder = T0Folder;
FtName = 'Ft.xlsx';
FtPath = fullfile(FtFolder, FtName);

F0aFolder = T0Folder;
F0aName = 'Fn_A.xlsx';
F0aPath = fullfile(F0aFolder, F0aName);

F0bFolder = T0Folder;
F0bName = 'FnQ_B.xlsx';
F0bPath = fullfile(F0bFolder, F0bName);

saveFolder = 'pont sens/article analyse donnees reelles/save';
saveName = ['TGV', num2str(N), 'save'];
if mode == 3
    saveName = [saveName, 'A'];
elseif mode == 4 || mode == 5
    saveName = [saveName, 'B'];
end
savePath = fullfile(saveFolder, saveName);


%% data loading

% t0
T0 = readtable(T0Path);
t0 = T0.Var1(N);

% ft, train freq estimates
Ft = readtable(FtPath);
ft = Ft.Var1(N);

% f0, Q, natural freq estimates (A) & Qmin
F0Qa = readtable(F0aPath);

% f0, Q, natural freq estimates (after) & Qmin
F0Qb = readtable(F0bPath);
f0b = F0Qb.f0;
Q0b = F0Qb.Q;

% acceleration
load(dataPath);
X = transpose(X);
t = X(1, :);
X = X(channels, :);
X = X - mean(X, 2);

% correction du temps
dt = diff(t);
dt = dt(abs(dt) < 0.25);
if any( abs(dt/mean(dt) - 1) > 1e-3)
    warning('pas de temps non constant');
end
dt = mean(dt);
t = dt * (0:length(t)-1);


%% premier affichage

if mode == 0
    figure;
    plt = plot(t, X);
    WaveletMenu('WaveletPlot', plt, 'fmin', 2, 'fmax', 20, 'Q', 10,...
        'RealShapePlot', @realShapePlotPont, 'MultiSignalMode', true, 'MultiSignalModeAbsValue', true,...
        'StopRidgeWhenIncreasing', true);
end


%% t0

if mode == 1
    figure;
    plt = plot(t, X);
end


%% t0 affichage

if mode == 1.5
    fig = figure;
    plt = plot(t, X);
    xlabel('Time [s]');
    ylabel('Acceleration [m/s²]');
    
    v_t = L_w * ft;
    Delta_t = (L_p + L_loco)/v_t;
    t1 = t0;
    t2 = t0 + (L_p + L_t)/v_t;
    
    xline(t1, '--', 'T_1', 'LabelHorizontalAlignment', 'left', 'LabelOrientation', 'horizontal');
    xline(t1+Delta_t, '--', 'T_1+\DeltaT', 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal');
    xline(t2-Delta_t, '--', {'', 'T_2-\DeltaT'}, 'LabelHorizontalAlignment', 'left', 'LabelOrientation', 'horizontal');
    xline(t2, '--', {'', 'T_2'}, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal');
    
    strLegend = [repmat('ch', 8, 1), num2str((1:8).')];
    legend(strLegend, 'FontSize', 8);
    
    fig.Position(3:4) = 8/10 * [560, 420];
    set(gca, 'YLim', max(abs(get(gca, 'YLim'))) * [-1 1]);
    set(gca, 'XLim', [46 55]);
end


%% ft

if mode == 2
    if ft == 0
        ft_estim = ft0;
    else
        ft_estim = ft;
    end
    disp(['ft_estim = ', num2str(ft_estim)]);
    v_t = L_w * ft_estim;
    Delta_t = (L_p + L_loco)/v_t;
    t1 = t0;
    t2 = t0 + (L_p + L_t)/v_t;
    
    figure;
    plt = plot(t, X);
    WaveletMenu('WaveletPlot', plt, 'fmin', 2, 'fmax', 7, 'Q', 7, 'XLim', [t1+Delta_t, t2-Delta_t],...
        'RealShapePlot', @realShapePlotPont, 'MultiSignalMode', true, 'MultiSignalModeAbsValue', true,...
        'StopRidgeWhenIncreasing', true);
end

%% extraction informations modales 1

if mode == 3
    v_t = L_w * ft;
    Delta_t = (L_p + L_loco)/v_t;
    t1 = t0;
    t2 = t0 + (L_p + L_t)/v_t;
    
    figure;
    plt = plot(t, X);
    
    shapes = nan(1, 8);
    shapesFt = nan(2, 8);
    ampl_graph = cell(1);
    freq_graph = cell(1);
    
    
    % mode 1
    disp(' ');
    disp('~~~~ mode 1 ~~~~');
    
    wvltFig = WaveletMenu('WaveletPlot', plt, 'fmin', 4.6, 'fmax', 7.5, 'Q', F0Qa.Qf1(N),...
        'XLim', [t1+Delta_t, t2-Delta_t],...
        'RealShapePlot', realShapePlotPont, 'MultiSignalMode', true, 'MultiSignalModeAbsValue', true,...
        'StopRidgeWhenIncreasing', false);
    
    % ampl-freq graph
    disp('ampl-freq graph');
    [ampl_graph{1}, freq_graph{1}] = PlotExtract();
    
    % mode shapes
    try
        shapes(1, 1) = input('shape: ');
    catch
        shapes(1, 1) = nan;
    end
    if ~isnan(shapes(1, 1))
        for kch = 2:8
            shapes(1, kch) = input('');
        end
    end
    % fin
    waitfor(wvltFig);
    
    
    % ft
    disp(' ');
    disp('~~~~ ft ~~~~');
    
    wvltFig = WaveletMenu('WaveletPlot', plt, 'fmin', 3, 'fmax', 5.5, 'Q', F0Qa.Qft(N),...
        'XLim', [t1+Delta_t, t2-Delta_t],...
        'RealShapePlot', realShapePlotPont, 'MultiSignalMode', true, 'MultiSignalModeAbsValue', true,...
        'StopRidgeWhenIncreasing', false);
    % "mode" shapes
    try
        shapesFt(1, 1) = input('shape: ');
    catch
        shapesFt(1, 1) = nan;
    end
    if ~isnan(shapesFt(1, 1))
        for kch = 2:8
            shapesFt(1, kch) = input('');
        end
    end
    % fin
    waitfor(wvltFig);
    
    
    % 2ft
    disp(' ');
    disp('~~~~ 2ft ~~~~');
    
    wvltFig = WaveletMenu('WaveletPlot', plt, 'fmin', 7.5, 'fmax', 9, 'Q', F0Qa.Qft2(N),...
        'XLim', [t1+Delta_t, t2-Delta_t],...
        'RealShapePlot', realShapePlotPont, 'MultiSignalMode', true, 'MultiSignalModeAbsValue', true,...
        'StopRidgeWhenIncreasing', false);
    % "mode" shapes
    try
        shapesFt(2, 1) = input('shape: ');
    catch
        shapesFt(2, 1) = nan;
    end
    if ~isnan(shapesFt(2, 1))
        for kch = 2:8
            shapesFt(2, kch) = input('');
        end
    end
    % fin
    waitfor(wvltFig);
    
    
    
    
    if saveResults
        save(savePath, 'shapes', 'shapesFt', 'ampl_graph', 'freq_graph');
        disp('results saved');
    end
    
end


%% extraction informations modales 2

if mode == 4
    v_t = L_w * ft;
    Delta_t = (L_p + L_loco)/v_t;
    t1 = t0;
    t2 = t0 + (L_p + L_t)/v_t;
    
    figure;
    plt = plot(t, X);
    
    Nmodes = length(f0b);
    shapes = nan(Nmodes, 8);
    ampl_graph = cell(1, Nmodes);
    freq_graph = cell(1, Nmodes);
    time_graph2 = cell(1, Nmodes);
    ampl_graph2 = cell(1, Nmodes);
    for k = 1:Nmodes
        disp(' ');
        disp(['~~~~ mode ', num2str(k), ' ~~~~']);
        
        fmin = f0b(k) - 2.5;
        fmax = f0b(k) + 4;
        if k == 2
            fmin = f0b(k) - 0.3;
        elseif k == 3
            fmax = f0b(k) +3;
        end
        wvltFig = WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q0b(k),...
            'XLim', [t2, min(t(end), t2+10)],...
            'RealShapePlot', realShapePlotPont, 'MultiSignalMode', true, 'MultiSignalModeAbsValue', true,...
            'StopRidgeWhenIncreasing', true);
        
        % time-ampl graph
        disp('time-ampl graph');
        [time_graph2{k}, ampl_graph2{k}] = PlotExtract();
        
        % ampl-freq graph
        disp('ampl-freq graph');
        [ampl_graph{k}, freq_graph{k}] = PlotExtract();
        
        % mode shapes
        try
            shapes(k, 1) = input('shape: ');
        catch
            shapes(k, 1) = nan;
        end
        if ~isnan(shapes(k, 1))
            for kch = 2:8
                shapes(k, kch) = input('');
            end
        end
        
        % fin
        waitfor(wvltFig);
    end
    
    if saveResults
        save(savePath, 'shapes', 'ampl_graph', 'freq_graph', 'time_graph2', 'ampl_graph2');
        disp('results saved');
    end
    
end















