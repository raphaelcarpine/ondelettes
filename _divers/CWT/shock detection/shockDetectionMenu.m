function shockDetectionMenu(t, X, signalChannels,...
    QMeanInit, MotherWaveletMeanInit, ctEdgeEffectsMeanInit, QSpectrumInit, MotherWaveletSpectrumInit,...
    freqLimMeanInit, nbFreqMeanInit, scaleFreqMeanInit, scaleMeanInit,...
    freqLimSpectrumInit, nbFreqSpectrumInit, scaleFreqSpectrumInit, scaleSpectrumInit,...
    multiChannelMeanInit, multiChannelSpectrumInit)
%SHOCKDETECTIONMENU Summary of this function goes here
%   Detailed explanation goes here

% parametres par defaut
meanFuncInit = @(x) abs(x).^2;
thresholdModeInit = 'mean'; % 'mean', 'absolute'
thresholdValueInit = 2;
maxDetectionMethodInit = 'local'; % 'local', 'global'  (detect local maxima or global maxima)
plotMeanInit = true;
plotSpectrumInit = true;
averageSpectrumInit = false;
% computeWholeWaveletInit = false; % à coder ?

NpcPCAInit = 2;
logScalePCAInit = false;
stdScalePCAInit = false;
plotScatterPCAInit = true;
plotPCPCAInit = true;
plotDistribPCAInit = true;

if nargin == 0 % test
    t = linspace(0, 40, 10000);
    X = sin(2*pi*6*(t-8)) .* exp(-0.02*2*pi*6*(t-8)) .* (t >= 8)...
        + 0.5 * sin(2*pi*6*(t-25)) .* exp(-0.02*2*pi*6*(t-25)) .* (t >= 25);
    X = [1; -2] * X;
    figure;
    plot(t, X);
    signalChannels = 3 + (1:size(X, 1));
    
    QMeanInit = 5;
    MotherWaveletMeanInit = 'morlet';
    ctEdgeEffectsMeanInit = [3 3];
    QSpectrumInit = 20;
    MotherWaveletSpectrumInit = 'cauchy';
    freqLimMeanInit = [1 10];
    nbFreqMeanInit = 300;
    scaleFreqMeanInit = 'log';
    scaleMeanInit = 'lin';
    freqLimSpectrumInit = [1 20];
    nbFreqSpectrumInit = 300;
    scaleFreqSpectrumInit = 'log';
    scaleSpectrumInit = 'lin';
    multiChannelMeanInit = true;
    multiChannelSpectrumInit = false;
end

%% input array size

% array size
if size(t, 1) == 1 && size(t, 2) == size(X, 2)
    % ok
elseif size(t, 1) == size(X, 1) && size(t, 2) == size(X, 2)
    dt = mean(diff(t(1, :)));
    for k_t = 2:size(t, 1)
        if max(abs(t(k_t, :) - t(1, :))) > dt * 1e-3
            error('different time arrays');
        end
    end
    t = t(1, :);
else
    error(['array size problem (', num2str(size(t)), ' & ', num2str(size(X)), ')']);
end

% time step
if any(abs(diff(t)/mean(diff(t)) - 1) > 1e-3)
    error('non-constant time step');
end

%% input window

HpansShock = [2, 5, 5, 9, 5, 2.5];
HpansSpectrum = [2, 5, 5, 7.5, 2.5];
HpansPCA = [9, 8.5, 2.5];

Hpans = [sum(HpansShock), sum(HpansSpectrum); 0, sum(HpansPCA)];

fig = figure('Name', 'Shock Detection Menu', 'numbertitle', 'off');
fig.Units = 'characters';
fig.Position(3) = 2*75;
fig.Position(2) = fig.Position(2) - max(sum(Hpans, 1)) + fig.Position(4);
fig.Position(4) = max(sum(Hpans, 1));
fig.MenuBar = 'none';
fig.Resize = false;


shockPan = uipanel('Parent', fig, 'Units', 'normalized', 'Title', 'Shock detection');
spectrumPan = uipanel('Parent', fig, 'Units', 'normalized', 'Title', 'Spectrums');
PCAPan = uipanel('Parent', fig, 'Units', 'normalized', 'Title', 'Principal component analysis');

margin = 0.01;
hpans = Hpans/max(sum(Hpans, 1)) - 2*margin;
for iz = 1:size(hpans, 1)
    for jz = 1:size(hpans, 2)
        if iz == 1
            zpans(iz, jz) = 1 - hpans(iz, jz) - margin;
        else
            zpans(iz, jz) = zpans(iz-1, jz) - hpans(iz, jz) - 2*margin;
        end
    end
end

shockPan.Position = [0.01, zpans(1, 1), 0.485, hpans(1, 1)];
spectrumPan.Position = [0.505, zpans(1, 2), 0.485, hpans(1, 2)];
PCAPan.Position = [0.505, zpans(2, 2), 0.485, hpans(2, 2)];

% shock panel

freqShockPan = uipanel('Parent',shockPan, 'Units', 'normalized', 'Title', 'Frequencies');
wvltShockPan = uipanel('Parent',shockPan, 'Units', 'normalized', 'Title', 'Wavelet');
shockShockPan = uipanel('Parent',shockPan, 'Units', 'normalized', 'Title', 'Shock detection');
plotShockPan = uipanel('Parent',shockPan, 'Units', 'normalized', 'Title', 'Plot');

margin = 0.02;
hpansShock = (1 - margin*(length(HpansShock)+1)) * HpansShock/sum(HpansShock);
zpansShock = margin * ones(size(hpansShock));
for iz = length(zpansShock)-1:-1:1
    zpansShock(iz) = zpansShock(iz+1) + hpansShock(iz+1) + margin;
end

freqShockPan.Position = [0.02, zpansShock(2), 0.96, hpansShock(2)];
wvltShockPan.Position = [0.02, zpansShock(3), 0.96, hpansShock(3)];
shockShockPan.Position = [0.02, zpansShock(4), 0.96, hpansShock(4)];
plotShockPan.Position = [0.02, zpansShock(5), 0.96, hpansShock(5)];

% spectrum panel

freqSpectrumPan = uipanel('Parent',spectrumPan, 'Units', 'normalized', 'Title', 'Frequencies');
wvltSpectrumPan = uipanel('Parent',spectrumPan, 'Units', 'normalized', 'Title', 'Wavelet');
plotSpectrumPan = uipanel('Parent',spectrumPan, 'Units', 'normalized', 'Title', 'Plot');

margin = 0.02;
hpansSpectrum = (1 - margin*(length(HpansSpectrum)+1)) * HpansSpectrum/sum(HpansSpectrum);
zpansSpectrum = margin * ones(size(hpansSpectrum));
for iz = length(zpansSpectrum)-1:-1:1
    zpansSpectrum(iz) = zpansSpectrum(iz+1) + hpansSpectrum(iz+1) + margin;
end

freqSpectrumPan.Position = [0.02, zpansSpectrum(2), 0.96, hpansSpectrum(2)];
wvltSpectrumPan.Position = [0.02, zpansSpectrum(3), 0.96, hpansSpectrum(3)];
plotSpectrumPan.Position = [0.02, zpansSpectrum(4), 0.96, hpansSpectrum(4)];

% PCA panel

paramPCAPan = uipanel('Parent', PCAPan, 'Units', 'normalized', 'Title', 'Parameters');
plotPCAPan = uipanel('Parent', PCAPan, 'Units', 'normalized', 'Title', 'Plot');

margin = 0.02;
hpansPCA = (1 - margin*(length(HpansPCA)+1)) * HpansPCA/sum(HpansPCA);
zpansPCA = margin * ones(size(hpansPCA));
for iz = length(zpansPCA)-1:-1:1
    zpansPCA(iz) = zpansPCA(iz+1) + hpansPCA(iz+1) + margin;
end

paramPCAPan.Position = [0.02, zpansPCA(1), 0.96, hpansPCA(1)];
plotPCAPan.Position = [0.02, zpansPCA(2), 0.96, hpansPCA(2)];

%% shock panel

% multi channel mode
multiChannelMeanInput = uicontrol('Parent', shockPan, 'Style', 'checkbox', 'String', 'single detection (multiple chanels)',...
    'Value', multiChannelMeanInit, 'Units', 'normalized', 'Position', [0.02, zpansShock(1), 0.97, hpansShock(1)]);

% freq panel
uicontrol('Parent', freqShockPan, 'Style', 'text', 'String', 'fmin:',...
    'Units', 'normalized', 'Position', [0., 0.15, 0.1, 0.6]);
fminMeanInput = uicontrol('Parent', freqShockPan, 'Style', 'edit', 'String', freqLimMeanInit(1),...
    'Units', 'normalized', 'Position', [0.1, 0.2, 0.1, 0.6]);
uicontrol('Parent', freqShockPan, 'Style', 'text', 'String', 'fmax:',...
    'Units', 'normalized', 'Position', [0.23, 0.15, 0.1, 0.6]);
fmaxMeanInput = uicontrol('Parent', freqShockPan, 'Style', 'edit', 'String', freqLimMeanInit(2),...
    'Units', 'normalized', 'Position', [0.33, 0.2, 0.1, 0.6]);
uicontrol('Parent', freqShockPan, 'Style', 'text', 'String', 'nb freqs:',...
    'Units', 'normalized', 'Position', [0.46, 0.15, 0.15, 0.6]);
nbFreqMeanInput = uicontrol('Parent', freqShockPan, 'Style', 'edit', 'String', nbFreqMeanInit,...
    'Units', 'normalized', 'Position', [0.61, 0.2, 0.1, 0.6]);
uicontrol('Parent', freqShockPan, 'Style', 'text', 'String', 'scale:',...
    'Units', 'normalized', 'Position', [0.74, 0.15, 0.1, 0.6]);
scaleFreqMeanValues = {'lin', 'log'};
scaleFreqMeanInput = uicontrol('Parent', freqShockPan, 'Style', 'popup', 'String', scaleFreqMeanValues,...
    'Value', find(strcmp(scaleFreqMeanValues, scaleFreqMeanInit)),...
    'Units', 'normalized', 'Position', [0.84, 0.27, 0.15, 0.6]);

% Q panel
uicontrol('Parent', wvltShockPan, 'Style', 'text', 'String', 'Q:',...
    'Units', 'normalized', 'Position', [0., 0.15, 0.1, 0.6]);
QMeanInput = uicontrol('Parent', wvltShockPan, 'Style', 'edit', 'String', QMeanInit,...
    'Units', 'normalized', 'Position', [0.1, 0.2, 0.1, 0.6]);
uicontrol('Parent', wvltShockPan, 'Style', 'text', 'String', 'mother wavelet:',...
    'Units', 'normalized', 'Position', [0.25, 0.15, 0.25, 0.6]);
MotherWaveletMeanValues = {'cauchy', 'morlet', 'harmonic', 'littlewood-paley', 'exponential'};
MotherWaveletMeanInput = uicontrol('Parent', wvltShockPan, 'Style', 'popup', 'String', MotherWaveletMeanValues,...
    'Value', find(strcmp(MotherWaveletMeanValues, MotherWaveletMeanInit)),...
    'Units', 'normalized', 'Position', [0.5, 0.27, 0.2, 0.6]);
uicontrol('Parent', wvltShockPan, 'Style', 'text', 'String', 'ct:',...
    'Units', 'normalized', 'Position', [0.76, 0.15, 0.05, 0.6]);
ctEdgeEffectsMeanInput = uicontrol('Parent', wvltShockPan, 'Style', 'edit', 'String', num2str(ctEdgeEffectsMeanInit),...
    'Units', 'normalized', 'Position', [0.81, 0.2, 0.15, 0.6]);

% shock panel
uicontrol('Parent', shockShockPan, 'Style', 'text', 'String', 'averaging function:',...
    'Units', 'characters', 'Position', [1, 3.5, 20, 1.5]);
meanFuncInput = uicontrol('Parent', shockShockPan, 'Style', 'edit', 'String', func2str(meanFuncInit),...
    'Units', 'characters', 'Position', [21, 3.5, 40, 1.8]);

uicontrol('Parent', shockShockPan, 'Style', 'text', 'String', 'threshold:',...
    'Units', 'characters', 'Position', [1, 0.8, 11, 1.5]);
thresholdValueInput = uicontrol('Parent', shockShockPan, 'Style', 'edit', 'String', thresholdValueInit,...
    'Units', 'characters', 'Position', [12, 0.8, 6, 1.8]);
thresholdModeValues = {'mean', 'absolute'};
thresholdModeInput = uicontrol('Parent', shockShockPan, 'Style', 'popup', 'String', thresholdModeValues,...
    'Value', find(strcmp(thresholdModeValues, thresholdModeInit)),...
    'Units', 'characters', 'Position', [19, 0.8, 13, 1.8]);

uicontrol('Parent', shockShockPan, 'Style', 'text', 'String', 'max detection:',...
    'Units', 'characters', 'Position', [40, 0.8, 15, 1.5]);
maxDetectionMethodValues = {'local', 'global'};
maxDetectionMethodInput = uicontrol('Parent', shockShockPan, 'Style', 'popup', 'String', maxDetectionMethodValues,...
    'Value', find(strcmp(maxDetectionMethodValues, maxDetectionMethodInit)),...
    'Units', 'characters', 'Position', [55, 0.8, 11, 1.8]);

% plot panel
plotMeanInput = uicontrol('Parent', plotShockPan, 'Style', 'checkbox', 'String', 'plot average',...
    'Value', plotMeanInit, 'Units', 'characters', 'Position', [1, 0.7, 20, 1.5]);
uicontrol('Parent', plotShockPan, 'Style', 'text', 'String', 'amplitude scale:',...
    'Units', 'characters', 'Position', [27, 0.5, 17, 1.5]);
scaleMeanValues = {'lin', 'log'};
scaleMeanInput = uicontrol('Parent', plotShockPan, 'Style', 'popup', 'String', scaleMeanValues,...
    'Value', find(strcmp(scaleMeanValues, scaleMeanInit)),...
    'Units', 'characters', 'Position', [44, 0.5, 8, 1.8]);

% ok button
computeShockInput = uicontrol('Parent', shockPan, 'Style', 'pushbutton', 'String', 'compute',...
    'Units', 'normalized', 'Position', [0.3, zpansShock(6), 0.4, hpansShock(6)]);


%% spectrum panel

% multi channel
multiChannelSpectrumInput = uicontrol('Parent', spectrumPan, 'Style', 'checkbox', 'String', 'average squared spectrums (multiple chanels)',...
    'Value', multiChannelSpectrumInit, 'Units', 'normalized', 'Position', [0.02, zpansSpectrum(1), 0.97, hpansSpectrum(1)]);

% interraction mulit channel mean
    function multiChannelMeanCallback(flag)
        if flag
            set(multiChannelSpectrumInput, 'Enable', 'on');
        else
            set(multiChannelSpectrumInput, 'Enable', 'off');
            set(multiChannelSpectrumInput, 'Value', false);
        end
    end

multiChannelMeanCallback(multiChannelMeanInit);
multiChannelMeanInput.Callback = @(hObject, ~) multiChannelMeanCallback(get(hObject, 'Value'));

% freq panel
uicontrol('Parent', freqSpectrumPan, 'Style', 'text', 'String', 'fmin:',...
    'Units', 'normalized', 'Position', [0., 0.15, 0.1, 0.6]);
fminSpectrumInput = uicontrol('Parent', freqSpectrumPan, 'Style', 'edit', 'String', freqLimSpectrumInit(1),...
    'Units', 'normalized', 'Position', [0.1, 0.2, 0.1, 0.6]);
uicontrol('Parent', freqSpectrumPan, 'Style', 'text', 'String', 'fmax:',...
    'Units', 'normalized', 'Position', [0.23, 0.15, 0.1, 0.6]);
fmaxSpectrumInput = uicontrol('Parent', freqSpectrumPan, 'Style', 'edit', 'String', freqLimSpectrumInit(2),...
    'Units', 'normalized', 'Position', [0.33, 0.2, 0.1, 0.6]);
uicontrol('Parent', freqSpectrumPan, 'Style', 'text', 'String', 'nb freqs:',...
    'Units', 'normalized', 'Position', [0.46, 0.15, 0.15, 0.6]);
nbFreqSpectrumInput = uicontrol('Parent', freqSpectrumPan, 'Style', 'edit', 'String', nbFreqSpectrumInit,...
    'Units', 'normalized', 'Position', [0.61, 0.2, 0.1, 0.6]);
uicontrol('Parent', freqSpectrumPan, 'Style', 'text', 'String', 'scale:',...
    'Units', 'normalized', 'Position', [0.74, 0.15, 0.1, 0.6]);
scaleFreqSpectrumValues = {'lin', 'log'};
scaleFreqSpectrumInput = uicontrol('Parent', freqSpectrumPan, 'Style', 'popup', 'String', scaleFreqSpectrumValues,...
    'Value', find(strcmp(scaleFreqSpectrumValues, scaleFreqSpectrumInit)),...
    'Units', 'normalized', 'Position', [0.84, 0.27, 0.15, 0.6]);

% Q panel
uicontrol('Parent', wvltSpectrumPan, 'Style', 'text', 'String', 'Q:',...
    'Units', 'normalized', 'Position', [0., 0.15, 0.1, 0.6]);
QSpectrumInput = uicontrol('Parent', wvltSpectrumPan, 'Style', 'edit', 'String', QSpectrumInit,...
    'Units', 'normalized', 'Position', [0.1, 0.2, 0.1, 0.6]);
uicontrol('Parent', wvltSpectrumPan, 'Style', 'text', 'String', 'mother wavelet:',...
    'Units', 'normalized', 'Position', [0.25, 0.15, 0.25, 0.6]);
MotherWaveletSpectrumValues = {'cauchy', 'morlet', 'harmonic', 'littlewood-paley', 'exponential'};
MotherWaveletSpectrumInput = uicontrol('Parent', wvltSpectrumPan, 'Style', 'popup', 'String', MotherWaveletSpectrumValues,...
    'Value', find(strcmp(MotherWaveletSpectrumValues, MotherWaveletSpectrumInit)),...
    'Units', 'normalized', 'Position', [0.5, 0.27, 0.2, 0.6]);

% plot panel
plotSpectrumInput = uicontrol('Parent', plotSpectrumPan, 'Style', 'checkbox', 'String', 'plot spectrums',...
    'Value', plotSpectrumInit, 'Units', 'characters', 'Position', [1, 2.5, 20, 1.5]);
uicontrol('Parent', plotSpectrumPan, 'Style', 'text', 'String', 'amplitude scale:',...
    'Units', 'characters', 'Position', [35, 2.5, 17, 1.5]);
scaleSpectrumValues = {'lin', 'log'};
scaleSpectrumInput = uicontrol('Parent', plotSpectrumPan, 'Style', 'popup', 'String', scaleSpectrumValues,...
    'Value', find(strcmp(scaleSpectrumValues, scaleSpectrumInit)),...
    'Units', 'characters', 'Position', [52, 2.5, 8, 1.8]);
averageSpectrumInput = uicontrol('Parent', plotSpectrumPan, 'Style', 'checkbox', 'String', 'plot average squared spectrum (multiple shocks)',...
    'Value', averageSpectrumInit, 'Units', 'characters', 'Position', [1, 0.5, 60, 1.5]);

% ok button
computeSpectrumInput = uicontrol('Parent', spectrumPan, 'Style', 'pushbutton', 'String', 'compute',...
    'Units', 'normalized', 'Position', [0.3, zpansSpectrum(5), 0.4, hpansSpectrum(5)]);


%% PCA panel

% param panel
uicontrol(paramPCAPan, 'Style', 'text', 'String', 'number of principal components:',...
    'Units', 'characters', 'Position', [1, 4.1, 33, 1.5], 'HorizontalAlignment', 'left');
NpcPCAInput = uicontrol(paramPCAPan, 'Style', 'edit', 'String', NpcPCAInit,...
    'Units', 'characters', 'Position', [34, 4.2, 7, 1.6]);
logScalePCAInput = uicontrol(paramPCAPan, 'Style', 'checkbox', 'String', 'log scale (amplitude)',...
    'Units', 'characters', 'Position', [1, 2.5, 40, 1.6], 'Value', logScalePCAInit);
stdScalePCAInput = uicontrol(paramPCAPan, 'Style', 'checkbox', 'String', 'scale by std dev',...
    'Units', 'characters', 'Position', [1, 0.7, 40, 1.6], 'Value', stdScalePCAInit);

% plot panel
plotScatterPCAInput = uicontrol(plotPCAPan, 'Style', 'checkbox', 'String', 'plot scatter',...
    'Units', 'characters', 'Position', [1, 3.7, 40, 1.6], 'Value', plotScatterPCAInit);
plotPCPCAInput = uicontrol(plotPCAPan, 'Style', 'checkbox', 'String', 'plot principal components',...
    'Units', 'characters', 'Position', [1, 2.2, 40, 1.6], 'Value', plotPCPCAInit);
plotDistribPCAInput = uicontrol(plotPCAPan, 'Style', 'checkbox', 'String', 'plot variance distribution',...
    'Units', 'characters', 'Position', [1, 0.7, 40, 1.6], 'Value', plotDistribPCAInit);

% compute button
computePCAInput = uicontrol(PCAPan, 'Style', 'pushbutton', 'String', 'compute',...
    'Units', 'normalized', 'Position', [0.3, zpansPCA(3), 0.4, hpansPCA(3)]);

% max Npc far scatter == 3
    function NpcPCAInputCallback(~,~)
        if str2double(get(NpcPCAInput, 'String')) > 3
            set(plotScatterPCAInput, 'Value', false);
            set(plotScatterPCAInput, 'Enable', 'off');
        elseif str2double(get(NpcPCAInput, 'String')) <= 3
            set(plotScatterPCAInput, 'Enable', 'on');
        end
    end

NpcPCAInput.Callback = @NpcPCAInputCallback;

%% shocks return function

computeShockInput.Callback = @(~, ~) computeShockReturn();

flagNewCWTShock = true;
flagNewShockDetection = true;

fminMean = nan;
fmaxMean = nan;
nbFreqMean = nan;
scaleFreqMean = '';
QMean = nan;
MotherWaveletMean = '';
ctEdgeEffectsMean = nan;
meanFunc = @() 0;

multiChannelMean = -1;
thresholdValue = nan;
thresholdMode = '';
maxDetectionMethod = '';

freqsMean = [];
meanWvltTot = nan(size(X));

meanWvlt = [];
shockIndexes = {};
thresholdAbsoluteValues = {};

    function computeShockReturn()
        % cursor
        set(fig, 'pointer', 'watch');
        drawnow;
        
        %%%% input
        
        % multi channel mode
        multiChannelMeanNew = get(multiChannelMeanInput, 'Value');
        
        % freq panel
        fminMeanNew = str2double(get(fminMeanInput, 'String'));
        fmaxMeanNew = str2double(get(fmaxMeanInput, 'String'));
        nbFreqMeanNew = str2double(get(nbFreqMeanInput, 'String'));
        scaleFreqMeanNew = scaleFreqMeanValues{get(scaleFreqMeanInput, 'Value')};
        
        % wvlt panel
        QMeanNew = str2double(get(QMeanInput, 'String'));
        MotherWaveletMeanNew =  MotherWaveletMeanValues{get(MotherWaveletMeanInput, 'Value')};
        ctEdgeEffectsMeanNew = str2num(get(ctEdgeEffectsMeanInput, 'String'));
        
        % shock panel
        meanFuncNew = str2func(get(meanFuncInput, 'String'));
        thresholdValueNew = str2double(get(thresholdValueInput, 'String'));
        thresholdModeNew = thresholdModeValues{get(thresholdModeInput, 'Value')};
        maxDetectionMethodNew = maxDetectionMethodValues{get(maxDetectionMethodInput, 'Value')};
        
        % plot panel
        plotMean = get(plotMeanInput, 'Value');
        scaleMean = scaleMeanValues{get(scaleMeanInput, 'Value')};
        
        
        %%%% check if new
        
        flagNewCWTShock = fminMeanNew ~= fminMean || fmaxMeanNew ~= fmaxMean || nbFreqMeanNew ~= nbFreqMean ||...
            ~strcmp(scaleFreqMeanNew, scaleFreqMean) || QMeanNew ~= QMean ||...
            ~strcmp(MotherWaveletMeanNew, MotherWaveletMean) || ~isequal(ctEdgeEffectsMeanNew, ctEdgeEffectsMean) ||...
            ~strcmp(func2str(meanFuncNew), func2str(meanFunc));
        
        flagNewShockDetection = multiChannelMean ~= multiChannelMeanNew || thresholdValue ~=thresholdValueNew ||...
            ~strcmp(thresholdMode, thresholdModeNew) || ~strcmp(maxDetectionMethod, maxDetectionMethodNew);
        
        flagNewShockDetection = flagNewShockDetection || flagNewCWTShock;
        
        fminMean = fminMeanNew;
        fmaxMean = fmaxMeanNew;
        nbFreqMean = nbFreqMeanNew;
        scaleFreqMean = scaleFreqMeanNew;
        QMean = QMeanNew;
        MotherWaveletMean = MotherWaveletMeanNew;
        ctEdgeEffectsMean = ctEdgeEffectsMeanNew;
        meanFunc = meanFuncNew;
        
        multiChannelMean = multiChannelMeanNew;
        thresholdValue = thresholdValueNew;
        thresholdMode = thresholdModeNew;
        maxDetectionMethod = maxDetectionMethodNew;
        
        
        %%%% compute CWT
        
        if flagNewCWTShock
            % frequencies array
            switch scaleFreqMean
                case 'lin'
                    freqsMean = linspace(fminMean, fmaxMean, nbFreqMean);
                case 'log'
                    freqsMean = logspace(log10(fminMean), log10(fmaxMean), nbFreqMean);
            end
            
            % computing meanWvltTot
            meanWvltTot = nan(size(X));
            for k_x = 1:size(X, 1)
                meanWvltTot(k_x, :) = WvltComp(t, X(k_x, :), freqsMean, QMean, 'MotherWavelet', MotherWaveletMean,...
                    'MeanOverFreqFunc', meanFunc);
            end
        end
        
        
        %%%% shock detection
        
        if multiChannelMean
            meanWvlt = mean(meanWvltTot, 1);
        else
            meanWvlt = meanWvltTot;
        end
        
        shockIndexes = cell(1, size(meanWvlt, 1));
        thresholdAbsoluteValues = cell(1, size(meanWvlt, 1));
        for k_ch = 1:size(meanWvlt, 1)
            figName = ['mean fcn: ', func2str(meanFunc), ' ; ',...
                num2str(fminMean), ' < f < ', num2str(fmaxMean), ' (', scaleFreqMean, ' scale) ; '];
            if multiChannelMean
                figName = [figName, 'all selected channels'];
            else
                figName = [figName, sprintf('channel %d', signalChannels(k_ch))];
            end
            [shockIndexes{k_ch}, thresholdAbsoluteValues{k_ch}] = ...
                maxDetection(t, meanWvlt(k_ch, :), freqsMean, QMean, MotherWaveletMean, ctEdgeEffectsMean,...
                thresholdMode, thresholdValue, maxDetectionMethod, plotMean, scaleMean, figName);
        end
        
        %%%% cursor
        set(fig, 'pointer', 'arrow');
        drawnow;
        
    end



%% spectrums return function

computeSpectrumInput.Callback = @(~, ~) computeSpectrumReturn();

flagNewCWTSpectrums = true;

fminSpectrum = nan;
fmaxSpectrum = nan;
nbFreqSpectrum = nan;
scaleFreqSpectrum = '';
QSpectrum = nan;
MotherWaveletSpectrum =  '';

spectrumsTot0 = {};
averageShockSpectrum0 = {};
averageSpectrum0 = {};
averageUnderThresholdSpectrum0 = {};
averageAboveThresholdSpectrum0 = {};

freqsSpectrum = [];
spectrumsTot = {};


    function computeSpectrumReturn()
        % cursor
        set(fig, 'pointer', 'watch');
        drawnow;
        
        %%%% compute shocks
        computeShockReturn();
        
        %%%% input
        
        % multi channel mode
        multiChannelSpectrum = get(multiChannelSpectrumInput, 'Value');
        % multi channel mode mean
        multiChannelMean = get(multiChannelMeanInput, 'Value');
        
        % freq panel
        fminSpectrumNew = str2double(get(fminSpectrumInput, 'String'));
        fmaxSpectrumNew = str2double(get(fmaxSpectrumInput, 'String'));
        nbFreqSpectrumNew = str2double(get(nbFreqSpectrumInput, 'String'));
        scaleFreqSpectrumNew = scaleFreqSpectrumValues{get(scaleFreqSpectrumInput, 'Value')};
        
        % wvlt panel
        QSpectrumNew = str2double(get(QSpectrumInput, 'String'));
        MotherWaveletSpectrumNew =  MotherWaveletSpectrumValues{get(MotherWaveletSpectrumInput, 'Value')};
        
        % plot panel
        plotShockSpectrums = get(plotSpectrumInput, 'Value');
        scaleSpectrum = scaleSpectrumValues{get(scaleSpectrumInput, 'Value')};
        
        % average spectrum mode
        plotAverageShockSpectrum = get(averageSpectrumInput, 'Value');
        
        
        %%%% check if new
        
        flagNewCWTSpectrums = fminSpectrum ~= fminSpectrumNew ||...
            fmaxSpectrum ~= fmaxSpectrumNew ||...
            nbFreqSpectrum ~= nbFreqSpectrumNew ||...
            ~strcmp(scaleFreqSpectrum, scaleFreqSpectrumNew) ||...
            QSpectrum ~= QSpectrumNew ||...
            ~strcmp(MotherWaveletSpectrum, MotherWaveletSpectrumNew);
        
        flagNewCWTSpectrums = flagNewCWTSpectrums || flagNewShockDetection;
        
        fminSpectrum = fminSpectrumNew;
        fmaxSpectrum = fmaxSpectrumNew;
        nbFreqSpectrum = nbFreqSpectrumNew;
        scaleFreqSpectrum = scaleFreqSpectrumNew;
        QSpectrum = QSpectrumNew;
        MotherWaveletSpectrum = MotherWaveletSpectrumNew;
        
        
        %%%% compute spectrum
            
        % frequencies array
        switch scaleFreqSpectrum
            case 'lin'
                freqsSpectrum = linspace(fminSpectrum, fmaxSpectrum, nbFreqSpectrum);
            case 'log'
                freqsSpectrum = logspace(log10(fminSpectrum), log10(fmaxSpectrum), nbFreqSpectrum);
        end
        
        % under & above threshold indexes
        UnderThresholdAverageIndexes = cell(1, size(meanWvlt, 1));
        AboveThresholdAverageIndexes = cell(1, size(meanWvlt, 1));
        for k_ch = 1:size(meanWvlt, 1)
            UnderThresholdAverageIndexes{k_ch} = meanWvlt(k_ch, :) < thresholdAbsoluteValues{k_ch};
            AboveThresholdAverageIndexes{k_ch} = meanWvlt(k_ch, :) >= thresholdAbsoluteValues{k_ch};
        end
        
        % already computed CWT
        if flagNewCWTSpectrums
            argsSavedCWT = {};
        else
            argsSavedCWT = {'spectrumsTot0', spectrumsTot0, 'averageShockSpectrum0', averageShockSpectrum0,...
                'averageSpectrum0', averageSpectrum0, 'averageUnderThresholdSpectrum0', averageUnderThresholdSpectrum0,...
                'averageAboveThresholdSpectrum0', averageAboveThresholdSpectrum0};
        end
        
        [spectrumsTot, spectrumsTot0, averageShockSpectrum0, averageSpectrum0,...
            averageUnderThresholdSpectrum0, averageAboveThresholdSpectrum0]...
            = getSpectrums(t, X, signalChannels, shockIndexes, freqsSpectrum, QSpectrum, MotherWaveletSpectrum,...
            scaleFreqSpectrum, scaleSpectrum,...
            multiChannelMean, multiChannelSpectrum, plotShockSpectrums,...
            'plotAverageShockSpectrum', plotAverageShockSpectrum,...
            'plotAverageSpectrum', plotAverageShockSpectrum,...
            'plotUnderThresholdAverage', plotAverageShockSpectrum,...
            'plotUnderThresholdAverageIndexes', UnderThresholdAverageIndexes,...
            'plotAboveThresholdAverage', plotAverageShockSpectrum,...
            'plotAboveThresholdAverageIndexes', AboveThresholdAverageIndexes,...
            argsSavedCWT{:});
        
        
        %%%% cursor
        set(fig, 'pointer', 'arrow');
        drawnow;
        
    end




%% PCA return function

computePCAInput.Callback = @(~, ~) computePCAReturn();

    function computePCAReturn()
        % cursor
        set(fig, 'pointer', 'watch');
        drawnow;
        
        %%%% compute spectrums
        computeSpectrumReturn();
        
        
        %%%% input
        
        % multi channel spectrum
        multiChannelSpectrum = get(multiChannelSpectrumInput, 'Value');
        % param panel
        NpcPCA = str2double(get(NpcPCAInput, 'String'));
        logScalePCA = get(logScalePCAInput, 'Value');
        stdScalePCA = get(stdScalePCAInput, 'Value');
        % plot panel
        plotScatterPCA = get(plotScatterPCAInput, 'Value');
        plotPCPCA = get(plotPCPCAInput, 'Value');
        plotDistribPCA = get(plotDistribPCAInput, 'Value');
        
        
        %%%% compute PCA
        
        for k_s = 1:length(spectrumsTot)
            if multiChannelSpectrum
                plotTitleSuffix = 'all selected channels';
            else
                plotTitleSuffix = sprintf('channel %u', signalChannels(k_s));
            end
            getPCA(spectrumsTot{k_s}, NpcPCA, logScalePCA, stdScalePCA, freqsSpectrum, scaleFreqSpectrum,...
                plotScatterPCA, plotPCPCA, plotDistribPCA, 'plotTitleSuffix', plotTitleSuffix, 'QSpectrum', QSpectrum);
        end
        
        
        %%%% cursor
        set(fig, 'pointer', 'arrow');
        drawnow;
        
    end

end

