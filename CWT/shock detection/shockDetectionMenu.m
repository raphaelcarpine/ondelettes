function shockDetectionMenu(t, X,...
    QMeanInit, MotherWaveletMeanInit, ctEdgeEffectsMeanInit, QSpectrumInit, MotherWaveletSpectrumInit,...
    freqLimMeanInit, nbFreqMeanInit, scaleFreqMeanInit, scaleMeanInit,...
    freqLimSpectrumInit, nbFreqSpectrumInit, scaleFreqSpectrumInit, scaleSpectrumInit,...
    multiSignalMeanInit, multiSignalSpectrumInit)
%SHOCKDETECTIONMENU Summary of this function goes here
%   Detailed explanation goes here

% parametres par defaut
meanFuncInit = @(x) abs(x).^2;
thresholdModeInit = 'mean'; % 'mean', 'absolute'
thresholdValueInit = 2;
maxDetectionMethodInit = 'local'; % 'local', 'global'  (detect local maxima or global maxima)
plotMeanInit = true;
plotSpectrumInit = true;
% computeWholeWaveletInit = false; % à coder ?


if nargin == 0 % test
    t = linspace(0, 40, 10000);
    X = sin(2*pi*6*(t-8)) .* exp(-0.02*2*pi*6*(t-8)) .* (t >= 8)...
        + 0.5 * sin(2*pi*6*(t-25)) .* exp(-0.02*2*pi*6*(t-25)) .* (t >= 25);
    X = [1; -2] * X;
    figure;
    plot(t, X);
    
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
    multiSignalMeanInit = true;
    multiSignalSpectrumInit = false;
end

%% input window

HpansShock = [5, 5, 9, 5, 2];
HpansSpectrum = [5, 5, 5, 2];

Hpans = [sum(HpansShock), sum(HpansSpectrum), 2];

fig = figure('Name', 'Shock Detection Menu', 'numbertitle', 'off');
fig.Units = 'characters';
fig.Position(3) = 75;
fig.Position(2) = fig.Position(2) - sum(Hpans) + fig.Position(4);
fig.Position(4) = sum(Hpans);
fig.MenuBar = 'none';


shockPan = uipanel('Parent',fig, 'Units', 'normalized', 'Title', 'Shock detection');
spectrumPan = uipanel('Parent',fig, 'Units', 'normalized', 'Title', 'Spectrums');

margin = 0.02;
hpans = (1 - margin*(length(Hpans)+1)) * Hpans/sum(Hpans);
zpans = margin * ones(size(hpans));
for iz = length(zpans)-1:-1:1
    zpans(iz) = zpans(iz+1) + hpans(iz+1) + margin;
end

shockPan.Position = [0.02, zpans(1), 0.96, hpans(1)];
spectrumPan.Position = [0.02, zpans(2), 0.96, hpans(2)];

% shock panel

freqShockPan = uipanel('Parent',shockPan, 'Units', 'normalized', 'Title', 'Frequencies');
wvltShockPan = uipanel('Parent',shockPan, 'Units', 'normalized', 'Title', 'Wavelet');
shockShockPan = uipanel('Parent',shockPan, 'Units', 'normalized', 'Title', 'Shock detection');
plotShockPan = uipanel('Parent',shockPan, 'Units', 'normalized', 'Title', 'Plot');

margin = 0.02;
hpans = (1 - margin*(length(HpansShock)+1)) * HpansShock/sum(HpansShock);
zpans = margin * ones(size(hpans));
for iz = length(zpans)-1:-1:1
    zpans(iz) = zpans(iz+1) + hpans(iz+1) + margin;
end

freqShockPan.Position = [0.02, zpans(1), 0.96, hpans(1)];
wvltShockPan.Position = [0.02, zpans(2), 0.96, hpans(2)];
shockShockPan.Position = [0.02, zpans(3), 0.96, hpans(3)];
plotShockPan.Position = [0.02, zpans(4), 0.96, hpans(4)];

% spectrum panel

freqSpectrumPan = uipanel('Parent',spectrumPan, 'Units', 'normalized', 'Title', 'Frequencies');
wvltSpectrumPan = uipanel('Parent',spectrumPan, 'Units', 'normalized', 'Title', 'Wavelet');
plotSpectrumPan = uipanel('Parent',spectrumPan, 'Units', 'normalized', 'Title', 'Plot');

margin = 0.02;
hpans = (1 - margin*(length(HpansSpectrum)+1)) * HpansSpectrum/sum(HpansSpectrum);
zpans = margin * ones(size(hpans));
for iz = length(zpans)-1:-1:1
    zpans(iz) = zpans(iz+1) + hpans(iz+1) + margin;
end

freqSpectrumPan.Position = [0.02, zpans(1), 0.96, hpans(1)];
wvltSpectrumPan.Position = [0.02, zpans(2), 0.96, hpans(2)];
plotSpectrumPan.Position = [0.02, zpans(3), 0.96, hpans(3)];

%% shock panel

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
uicontrol('Parent', plotShockPan, 'Style', 'text', 'String', 'amplitude sclae:',...
    'Units', 'characters', 'Position', [27, 0.5, 17, 1.5]);
scaleMeanValues = {'lin', 'log'};
scaleMeanInput = uicontrol('Parent', plotShockPan, 'Style', 'popup', 'String', scaleMeanValues,...
    'Value', find(strcmp(scaleMeanValues, scaleMeanInit)),...
    'Units', 'characters', 'Position', [44, 0.5, 8, 1.8]);

% multi signal mode
multiSignalMeanInput = uicontrol('Parent', shockPan, 'Style', 'checkbox', 'String', 'single detection (multiple chanels)',...
    'Value', multiSignalMeanInit, 'Units', 'characters', 'Position', [1, 0.5, 50, 1.5]);


%% spectrum panel

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
    'Value', plotSpectrumInit, 'Units', 'characters', 'Position', [1, 0.7, 20, 1.5]);
uicontrol('Parent', plotSpectrumPan, 'Style', 'text', 'String', 'amplitude sclae:',...
    'Units', 'characters', 'Position', [27, 0.5, 17, 1.5]);
scaleSpectrumValues = {'lin', 'log'};
scaleSpectrumInput = uicontrol('Parent', plotSpectrumPan, 'Style', 'popup', 'String', scaleSpectrumValues,...
    'Value', find(strcmp(scaleSpectrumValues, scaleSpectrumInit)),...
    'Units', 'characters', 'Position', [44, 0.5, 8, 1.8]);

% multi signal
multiSignalSpectrumInput = uicontrol('Parent', spectrumPan, 'Style', 'checkbox', 'String', 'sum of squared spectrums (multiple chanels)',...
    'Value', multiSignalSpectrumInit, 'Units', 'characters', 'Position', [1, 0.5, 50, 1.5]);

% interraction mulit signal mean
    function multiSignalMeanCallback(flag)
        if flag
            set(multiSignalSpectrumInput, 'Enable', 'on');
        else
            set(multiSignalSpectrumInput, 'Enable', 'off');
            set(multiSignalSpectrumInput, 'Value', false);
        end
    end

multiSignalMeanCallback(multiSignalMeanInit);
multiSignalMeanInput.Callback = @(hObject, ~) multiSignalMeanCallback(get(hObject, 'Value'));


%% ok button

okInput = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'OK',...
    'Value', find(strcmp(scaleSpectrumValues, scaleSpectrumInit)),...
    'Units', 'characters', 'Position', [30, 0.5, 15, 2.5]);

%% return func

okInput.Callback = @(~, ~) okReturn();

    function okReturn()
        %% input
        
        % freq panel
        fminMean = str2double(get(fminMeanInput, 'String'));
        fmaxMean = str2double(get(fmaxMeanInput, 'String'));
        nbFreqMean = str2double(get(nbFreqMeanInput, 'String'));
        scaleFreqMean = scaleFreqMeanValues{get(scaleFreqMeanInput, 'Value')};
        
        % wvlt panel
        QMean = str2double(get(QMeanInput, 'String'));
        MotherWaveletMean =  MotherWaveletMeanValues{get(MotherWaveletMeanInput, 'Value')};
        ctEdgeEffectsMean = str2num(get(ctEdgeEffectsMeanInput, 'String'));
        
        % shock panel
        meanFunc = str2func(get(meanFuncInput, 'String'));
        thresholdValue = str2double(get(thresholdValueInput, 'String'));
        thresholdMode = thresholdModeValues{get(thresholdModeInput, 'Value')};
        maxDetectionMethod = maxDetectionMethodValues{get(maxDetectionMethodInput, 'Value')};
        
        % plot panel
        plotMean = get(plotMeanInput, 'Value');
        scaleMean = scaleMeanValues{get(scaleMeanInput, 'Value')};
        
        % multi signal mode
        multiSignalMean = get(multiSignalMeanInput, 'Value');
        
        % freq panel
        fminSpectrum = str2double(get(fminSpectrumInput, 'String'));
        fmaxSpectrum = str2double(get(fmaxSpectrumInput, 'String'));
        nbFreqSpectrum = str2double(get(nbFreqSpectrumInput, 'String'));
        scaleFreqSpectrum = scaleFreqSpectrumValues{get(scaleFreqSpectrumInput, 'Value')};
        
        % wvlt panel
        QSpectrum = str2double(get(QSpectrumInput, 'String'));
        MotherWaveletSpectrum =  MotherWaveletSpectrumValues{get(MotherWaveletSpectrumInput, 'Value')};
        
        % plot panel
        plotSpectrum = get(plotSpectrumInput, 'Value');
        scaleSpectrum = scaleSpectrumValues{get(scaleSpectrumInput, 'Value')};
        
        % multi signal mode
        multiSignalSpectrum = get(multiSignalSpectrumInput, 'Value');
        
        
        %%
        multiSignalSpectrum = multiSignalSpectrum & multiSignalMean;
        
        
        %%
        
        switch scaleFreqMean
            case 'lin'
                freqsMean = linspace(fminMean, fmaxMean, nbFreqMean);
            case 'log'
                freqsMean = logspace(log10(fminMean), log10(fmaxMean), nbFreqMean);
        end
        
        switch scaleFreqSpectrum
            case 'lin'
                freqsSpectrum = linspace(fminSpectrum, fmaxSpectrum, nbFreqSpectrum);
            case 'log'
                freqsSpectrum = logspace(log10(fminSpectrum), log10(fmaxSpectrum), nbFreqSpectrum);
        end
        
        shockDetection(t, X, freqsMean, freqsSpectrum, QMean, QSpectrum, MotherWaveletMean, MotherWaveletSpectrum,...
            ctEdgeEffectsMean, meanFunc, thresholdMode, thresholdValue, maxDetectionMethod,...
            multiSignalMean, multiSignalSpectrum, 'plotMean', plotMean, 'plotSpectrum', plotSpectrum,...
            'meanScale', scaleMean, 'spectrumScale', scaleSpectrum);
        
        %%
        
%         delete(fig);
        
    end

end

