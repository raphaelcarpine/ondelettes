function fig = WaveletMenu(varargin)
%WaveletMenu Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

defaultX = nan;
defaultY = nan;
defaultFmin = 1;
defaultFmax = 10;
defaultNbFreq = 300;
defaultParent = 0;
defaultWaveletPlot = 0;
defaultQ = 10;
defaultMaxRidges = 1;
defaultMaxParallelRidges = inf;
defaultMaxSlopeRidge = 1;
defaultMultipleAxesDisplay = false;
defaultRidgeMinModu = 0;
defaultStopRidgeWhenIncreasing = false;
defaultCtEdgeEffects = 3;
defaultCf = 5;
defaultZeroPaddingFourier = 0;
defaultFourierAveraging = false;
defaultFourierAveragingNb = 10;
defaultFourierWindow = 'rectangular';
defaultFourierWindowParams = [];
defaultMultiSignalMode = false;
defaultMultiSignalModeAbsValue = false;
defaultRdtMode = false;
defaultAutocorrelationMode = false;
defaultAutocorrelationBias = 'biased'; %'biased', unbiased';
defaultAutocorrelationSVDMode = false;
defaultAutocorrelationFourierSVDMode = false;
defaultMaxLagCorr = nan;
defaultAutocorrelationNsvd = 1;
defaultWvltScale = 'log';
defaultFourierScale = 'lin';
defaultFrequencyScale = 'lin';
defaultXLim = nan;
defaultShowXLim = true;
defaultWvltAxesTitle = '';
defaultComplexShapePlot = @complexShapePlot1;
defaultRealShapePlot = @realShapePlot1;
defaultMotherWavelet = 'cauchy';
defaultXLimRidge = [-inf inf];
defaultctRidge = 0;
defaultRemoveMean = false;
defaultSignalUnit = 'm/s²';
defaultSquaredSignalUnit = 'm^2/s^4';
defaultCheckForUpdates = true;



checkParent = @(f) isa(f, 'matlab.ui.Figure') || isa(f, 'matlab.ui.container.Panel')...
    || isa(f, 'matlab.ui.container.Tab') || isa(f, 'matlab.ui.container.ButtonGroup');

addParameter(p,'WaveletPlot', defaultWaveletPlot); %si les données viennent d'une courbe directement (ou plusieurs)
addOptional(p, 'X', defaultX);
addOptional(p, 'Y', defaultY);
addParameter(p, 'fmin', defaultFmin);
addParameter(p, 'fmax', defaultFmax);
addParameter(p, 'NbFreq', defaultNbFreq);
addParameter(p,'Parent', defaultParent, checkParent);
addParameter(p,'Q', defaultQ);
addParameter(p,'MaxRidges', defaultMaxRidges);
addParameter(p,'MaxParallelRidges', defaultMaxParallelRidges);
addParameter(p,'MaxSlopeRidge', defaultMaxSlopeRidge);
addParameter(p,'MultipleAxesDisplay', defaultMultipleAxesDisplay);
addParameter(p,'RidgeMinModu', defaultRidgeMinModu);
addParameter(p,'StopRidgeWhenIncreasing', defaultStopRidgeWhenIncreasing);
addParameter(p,'CtEdgeEffects', defaultCtEdgeEffects);
addParameter(p,'Cf', defaultCf);
addParameter(p,'ZeroPaddingFourier', defaultZeroPaddingFourier);
addParameter(p,'FourierAveraging', defaultFourierAveraging);
addParameter(p,'FourierAveragingNb', defaultFourierAveragingNb);
addParameter(p,'FourierWindow', defaultFourierWindow);
addParameter(p,'FourierWindowParams', defaultFourierWindowParams);
addParameter(p,'MultiSignalMode', defaultMultiSignalMode);
addParameter(p,'MultiSignalModeAbsValue', defaultMultiSignalModeAbsValue);
addParameter(p,'RdtMode', defaultRdtMode);
addParameter(p,'AutocorrelationMode', defaultAutocorrelationMode);
addParameter(p,'AutocorrelationBias', defaultAutocorrelationBias);
addParameter(p,'AutocorrelationSVDMode', defaultAutocorrelationSVDMode);
addParameter(p,'AutocorrelationFourierSVDMode', defaultAutocorrelationFourierSVDMode);
addParameter(p,'AutocorrelationMaxLag', defaultMaxLagCorr);
addParameter(p,'AutocorrelationNsvd', defaultAutocorrelationNsvd);
addParameter(p,'WvltScale', defaultWvltScale);
addParameter(p,'FourierScale', defaultFourierScale);
addParameter(p,'FrequencyScale', defaultFrequencyScale);
addParameter(p,'XLim', defaultXLim);
addParameter(p,'ShowXLim', defaultShowXLim);
addParameter(p, 'WvltAxesTitle', defaultWvltAxesTitle);
addParameter(p, 'ComplexShapePlot', defaultComplexShapePlot);
addParameter(p, 'RealShapePlot', defaultRealShapePlot);
addParameter(p, 'MotherWavelet', defaultMotherWavelet);
addParameter(p, 'XLimRidge', defaultXLimRidge);
addParameter(p, 'ctRidge', defaultctRidge);
addParameter(p, 'RemoveMean', defaultRemoveMean);
addParameter(p, 'SignalUnits', defaultSignalUnit);
addParameter(p, 'SquaredSignalUnits', defaultSquaredSignalUnit);
addParameter(p, 'CheckForUpdates', defaultCheckForUpdates);

parse(p, varargin{:})

fig = p.Results.Parent;
if fig == 0
    fig = figure('Name', 'Wavelet Menu', 'numbertitle', 'off');
    fig.Units = 'characters';
    fig.Position(3) = 110;
    fig.Position(4) = 27;
    fig.MenuBar = 'none';
    fig.Resize = false;
%     fig.ToolBar = 'none';
end

    function M = cellmat2mat(M)
        if isa(M, 'cell')
            for iM = 1:length(M)
                M{iM} = M{iM}.';
            end
            M = [M{:}].';
        end
    end

    function closeParent(parent)
        try
            delete(fig);
        catch
        end
        delete(parent);
    end

WaveletPlot = p.Results.WaveletPlot;
if size(WaveletPlot, 1) > size(WaveletPlot, 2)
    WaveletPlot = transpose(WaveletPlot);
end

if WaveletPlot == 0
    getX = @() p.Results.X;
    getY = @() p.Results.Y;
    plotAxes = gobjects(1, 0);
    plotAxesName = '';
else
    getX = @() cellmat2mat (get(WaveletPlot, 'XData'));
    getY = @() cellmat2mat (get(WaveletPlot, 'YData'));
    plotAxes = cellmat2mat (get(WaveletPlot, 'Parent'));
    plotAxesName = ['fig', num2str(get(get(plotAxes(1), 'Parent'), 'number'))];
    set(fig, 'Name', [get(fig, 'Name'), ' (', plotAxesName, ')']);
    
    for waveplt = WaveletPlot
        parent = waveplt;
        while ~isa(parent, 'matlab.ui.Figure')
            parent = get(parent, 'Parent');
        end
        parent.CloseRequestFcn = @(~,~) closeParent(parent);
    end
end
nbPlotsTot = size(getY(), 1);
nbPlots = nbPlotsTot;
signalChannels = 1:nbPlots;


fmin = p.Results.fmin;
fmax = p.Results.fmax;
NbFreq = p.Results.NbFreq;
Q = p.Results.Q;
MaxRidges = p.Results.MaxRidges;
MaxParallelRidges = p.Results.MaxParallelRidges;
MaxSlopeRidge = p.Results.MaxSlopeRidge;
multipleAxesDisplay = p.Results.MultipleAxesDisplay;
RidgeMinModu = p.Results.RidgeMinModu;
StopRidgeWhenIncreasing = p.Results.StopRidgeWhenIncreasing;
ctEdgeEffects = p.Results.CtEdgeEffects;
cf = p.Results.Cf;
ZeroPaddingFourier = p.Results.ZeroPaddingFourier;
FourierAveraging = p.Results.FourierAveraging;
FourierAveragingNb = p.Results.FourierAveragingNb;
FourierWindow = p.Results.FourierWindow;
FourierWindowParams = p.Results.FourierWindowParams;

multiSignalMode = p.Results.MultiSignalMode;
multiSignalModeAbsValue = p.Results.MultiSignalModeAbsValue;
rdtMode = p.Results.RdtMode;
autocorrelationMode = p.Results.AutocorrelationMode;
autocorrelationBias = p.Results.AutocorrelationBias;
autocorrelationSVDMode = p.Results.AutocorrelationSVDMode;
autocorrelationFourierSVDMode = p.Results.AutocorrelationFourierSVDMode;
tRy = nan; Ry = nan; SVry = nan; SVvectry = nan;
resetCrossCorr();

autocorrelationNsvd = p.Results.AutocorrelationNsvd;
maxLagCorr = p.Results.AutocorrelationMaxLag;
WvltScale = p.Results.WvltScale;
FourierScale = p.Results.FourierScale;
FrequencyScale = p.Results.FrequencyScale;
XLim = p.Results.XLim;
XLimRidge = p.Results.XLimRidge;
ShowXLim = p.Results.ShowXLim;
wvltAxesTitle = p.Results.WvltAxesTitle;
ComplexShapePlot = p.Results.ComplexShapePlot;
RealShapePlot = p.Results.RealShapePlot;
MotherWavelet = p.Results.MotherWavelet;
ctRidge = p.Results.ctRidge;
removeMean = p.Results.RemoveMean;
signalUnit = p.Results.SignalUnits;
squaredSignalUnit = p.Results.SquaredSignalUnits;
CheckForUpdates = p.Results.CheckForUpdates;

x0 = getX();
if isnan(XLim)
    XLim = [x0(1), x0(end)];
end

if isnan(maxLagCorr)
    maxLagCorr = XLim(2)-XLim(1);
end

x = nan; y = nan;

%% reglage affichage subpolt/simple plot

if multipleAxesDisplay
    subplot0 = @(i,j,k,ax) subplot(i,j,k,ax);
else
    subplot0 = @(i,j,k,ax) ax;
end

    function setMultipleAxesDisplay(value)
        multipleAxesDisplay = value;
        if multipleAxesDisplay
            subplot0 = @(i,j,k,ax) subplot(i,j,k,ax);
        else
            subplot0 = @(i,j,k,ax) ax;
        end
    end

%% affichage XLim

XLimDisplayObj = XLimDisplay(XLim, XLimRidge, ShowXLim, ShowXLim, plotAxes, true(size(plotAxes)));

%% random decrement technique

Xrdt = [];
Yrdt = [];


%% bouton ondelettes et panneaux param et sorties

%ondelette
waveletPan = uipanel('Parent', fig, 'Units', 'normalized');
waveletPan.Position = [0 0.15 1 0.85];

buttonWavelet = uicontrol('Parent',waveletPan, 'Units', 'normalized','Style','pushbutton',...
    'String', 'wavelet');
paramPan = uipanel('Parent',waveletPan, 'Units', 'normalized', 'Title', 'parameters');
plotPan = uipanel('Parent',waveletPan, 'Units', 'normalized', 'Title', 'plots');
shapesPan = uipanel('Parent',waveletPan, 'Units', 'normalized', 'Title', 'mode shapes');

margin = 0.005;
vertProp = 0.15; % vertical proporitons
horProp = [0.25, 0.625]; % horizontal proportions
buttonWavelet.Position = [margin, margin, 1-2*margin, vertProp - 1.5*margin];
paramPan.Position = [margin, vertProp+0.5*margin, horProp(1)-1.5*margin, 1-vertProp-1.5*margin];
plotPan.Position = [horProp(1)+0.5*margin, vertProp+0.5*margin, horProp(2)-horProp(1)-margin, 1-vertProp-1.5*margin];
shapesPan.Position = [horProp(2)+0.5*margin, vertProp+0.5*margin, 1-horProp(2)-1.5*margin, 1-vertProp-1.5*margin];

%autres transformees
transformPan = uipanel('Parent',fig, 'Units', 'normalized');
transformPan.Position = [0 0 1 0.15];

buttonFourier = uicontrol('Parent',transformPan, 'Units', 'normalized','Style','pushbutton',...
    'String', 'fourier');
buttonHilbert = uicontrol('Parent',transformPan, 'Units', 'normalized','Style','pushbutton',...
    'String', 'hilbert');

buttonFourier.Position = [0.02, 0.01, 0.47, 0.98];
buttonHilbert.Position = [0.51, 0.01, 0.47, 0.98];


%% paramètres

strfmin = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'fmin = ');
editfmin = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(fmin));

strfmax = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'fmax = ');
editfmax = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(fmax));

strNbFreq = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'NbFreq = ');
editNbFreq = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(NbFreq));

strQ = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'Q = ');
editQ = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(Q));

strmaxR = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'max ridges');
editmaxR = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(MaxRidges));

strPR = uicontrol('Parent', paramPan, 'Units', 'normalized','Style','text',...
    'String', 'max parallel ridges');
editPR = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(MaxParallelRidges));

strSL = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'max ridge slope');
editSL = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(MaxSlopeRidge));

Strs = [strfmin, strfmax, strNbFreq, strQ, strmaxR, strPR, strSL];
Edits =[editfmin, editfmax, editNbFreq, editQ, editmaxR, editPR, editSL];
n = length(Strs);
for k=1:n
    Strs(k).Position = [0.01, 0.01+(n-k)/n, 0.48, 1/n-0.02];
    Edits(k).Position = [0.51, 0.01+(n-k)/n, 0.48, 1/n-0.02];
end



%% sorties

checkboxGeneral = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'general plot', 'Value', false);

checkboxModule = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'module plot', 'Value', false);
WvltScaleButton = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

checkboxPhase = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'phase plot', 'Value', false);

checkboxTimeAmplPlot = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'time, ampl on plot', 'Value', false);
timeAmplPlots = [];
deleteButton = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','pushbutton',...
    'String', 'delete');


Checkboxs1 = [checkboxGeneral, checkboxModule, checkboxPhase, checkboxTimeAmplPlot];


if isempty(plotAxes)
    set(checkboxTimeAmplPlot, 'Enable', 'off');
    set(deleteButton, 'Enable', 'off');
end





checkboxTimeAmpl = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'time, ampl', 'Value', false);
xscaleTimeAmpl = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleTimeAmpl = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

checkboxTimeFreq = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'time, freq', 'Value', false);
xscaleTimeFreq = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleTimeFreq = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

checkboxTimeDamp = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'time, damping', 'Value', false);
xscaleTimeDamp = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleTimeDamp = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

checkboxAmplFreq = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'ampl, freq', 'Value', false);
xscaleAmplFreq = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleAmplFreq = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

checkboxAmplDamp = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'ampl, damping', 'Value', false);
xscaleAmplDamp = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleAmplDamp = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

checkboxTimePhase = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'time, phase', 'Value', false);
xscaleTimePhase = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleTimePhase = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

Checkboxs2 = [checkboxTimeAmpl, checkboxTimeFreq, checkboxTimeDamp, checkboxAmplFreq, checkboxAmplDamp,...
    checkboxTimePhase];
XScales = [xscaleTimeAmpl, xscaleTimeFreq, xscaleTimeDamp, xscaleAmplFreq, xscaleAmplDamp, xscaleTimePhase];
YScales = [yscaleTimeAmpl, yscaleTimeFreq, yscaleTimeDamp, yscaleAmplFreq, yscaleAmplDamp, yscaleTimePhase];
n1 = length(Checkboxs1);
n2 = length(Checkboxs2);
n = n1+n2;
for k=1:n1
    Checkboxs1(k).Position = [0.01, 0.01+(n-k)/n, 0.48, 1/n-0.02];
end
Checkboxs1(end).Position(3) = 0.7;

for k=(n1+1:n1+n2)
    Checkboxs2(k-n1).Position = [0.01, 0.01+(n-k)/n, 0.48, 1/n-0.02];
end
WvltScaleButton.Position = [0.76, 0.01+(n-2)/n, 0.23, 1/n-0.02];
for k=2:n1
    deleteButton.Position = [0.76, 0.01+(n-k)/n, 0.23, 1/n-0.02];
end
for k=(n1+1:n1+n2)
    XScales(k-n1).Position = [0.51, 0.01+(n-k)/n, 0.23, 1/n-0.02];
    YScales(k-n1).Position = [0.76, 0.01+(n-k)/n, 0.23, 1/n-0.02];
end

scalesNames = {'linear', 'log'};
for k=(n1+1:n1+n2)
    XScales(k-n1).Callback = @(~,~) set(XScales(k-n1), 'String', scalesNames{XScales(k-n1).Value+1});
    YScales(k-n1).Callback = @(~,~) set(YScales(k-n1), 'String', scalesNames{YScales(k-n1).Value+1});
end

LogScales = [yscaleTimeAmpl, xscaleAmplFreq, xscaleAmplDamp];
for scale = LogScales
    set(scale, 'Value', 1);
    set(scale, 'String', scalesNames{scale.Value+1});
end

switch WvltScale
    case 'lin'
    case 'log'
        WvltScaleButton.Value = true;
        WvltScaleButton.String = 'log';
end

    function WvltScaleButtonCallback(~,~)
        if WvltScaleButton.Value
            WvltScale = 'log';
            WvltScaleButton.String = 'log';
        else
            WvltScale = 'lin';
            WvltScaleButton.String = 'linear';
        end
    end
WvltScaleButton.Callback = @WvltScaleButtonCallback;



%% sorties mode shapes

strShapesInstantaneous = uicontrol('Parent', shapesPan, 'Units', 'normalized', 'Style', 'text',...
    'String', 'instantaneous:', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

strxAxisShapes = uicontrol('Parent', shapesPan, 'Units', 'normalized', 'Style', 'text',...
    'String', 'x-axis:', 'HorizontalAlignment', 'left');
xAxisShapesnames = {'time', 'ampl', 'log(ampl)', 'freq'};
xAxisShapes = uicontrol('Parent', shapesPan, 'Units', 'normalized','Style', 'popupmenu',...
    'String', xAxisShapesnames);



stryAxisShapes = uicontrol('Parent', shapesPan, 'Units', 'normalized', 'Style', 'text',...
    'String', 'y-axis:', 'HorizontalAlignment', 'left');

checkboxRealShapes = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'real part', 'Value', false);
checkboxImagShapes = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'imag. part', 'Value', false);
checkboxModuleShapes = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'module', 'Value', false);
checkboxPhaseShapes = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'angle', 'Value', false);
checkboxAmplShapes = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'amplitude', 'Value', false);
Checkboxs3 = [checkboxRealShapes, checkboxImagShapes, checkboxModuleShapes, checkboxPhaseShapes, checkboxAmplShapes];



strShapesMean = uicontrol('Parent', shapesPan, 'Units', 'normalized', 'Style', 'text',...
    'String', 'average:', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
weightOptionShapesMean = uicontrol('Parent', shapesPan, 'Units', 'normalized', 'Style','togglebutton',...
    'String', 'weighted (ampl)', 'Value', false);

checkboxRealShapesMean = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'plot real', 'Value', false);
checkboxComplexShapesMean = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'plot complex', 'Value', false);
checkboxDispShapesMean = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'disp complex', 'Value', false);
checkboxDispFreqsMean = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'damped freq.', 'Value', false);
checkboxAmplRegMean = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'ampl. regression', 'Value', false);
Checkboxs4 = [checkboxRealShapesMean, checkboxComplexShapesMean, checkboxDispShapesMean,...
    checkboxDispFreqsMean, checkboxAmplRegMean];



n3 = ceil(length(Checkboxs3)/2);
n4 = ceil(length(Checkboxs4)/2);
n = n3 + n4 + 5.5;

strShapesInstantaneous.Position = [0.05, (n-1.2)/n, 0.94, 1/n-0.02];
strxAxisShapes.Position = [0.05, (n-2.2)/n, 0.94, 1/n-0.02];
xAxisShapes.Position = [0.01, (n-3)/n+0.01, 0.48, 1/n-0.02];

stryAxisShapes.Position = [0.05, (n-4.2)/n, 0.94, 1/n-0.02];
for k=1:length(Checkboxs3)
    Checkboxs3(k).Position = [0.5*mod(k-1, 2)+0.01, 0.01+(n-floor((k-1)/2)-5)/n, 0.48, 1/n-0.02];
end

strShapesMean.Position = [0.05, (n-5.5-n3)/n, 0.35, 1/n-0.02];
weightOptionShapesMean.Position = [0.37, (n-5.5-n3)/n, 0.52, 1/n-0.01];
for k=1:length(Checkboxs4)
    Checkboxs4(k).Position = [0.5*mod(k-1, 2)+0.01, 0.01+(n-floor((k-1)/2)-6.5-n3)/n, 0.48, 1/n-0.02];
end



%% menu (autres paramètres)

fig.MenuBar = 'none';

paramMenu = uimenu(fig,'Text','Options');

% mother wavelet
MotherWaveletNames = {'cauchy', 'morlet', 'harmonic', 'littlewood-paley', 'exponential'};

motherWaveletMenu = uimenu(paramMenu, 'Text','Mother wavelet');
motherWaveletMenuChoices(1) = uimenu(motherWaveletMenu, 'Text', 'Cauchy');
motherWaveletMenuChoices(2) = uimenu(motherWaveletMenu, 'Text', 'Morlet');
motherWaveletMenuChoices(3) = uimenu(motherWaveletMenu, 'Text', 'Harmonic');
motherWaveletMenuChoices(4) = uimenu(motherWaveletMenu, 'Text', 'Littlewood-Paley');
motherWaveletMenuChoices(5) = uimenu(motherWaveletMenu, 'Text', 'Exponential');
set(motherWaveletMenuChoices(find(strcmp(MotherWaveletNames, MotherWavelet))), 'Checked' ,'on');

    function selectMotherWaveletMenu(kchoice)
        for kchoices = 1:length(motherWaveletMenuChoices)
            set(motherWaveletMenuChoices(kchoices), 'Checked', 'off');
        end
        set(motherWaveletMenuChoices(kchoice), 'Checked', 'on');
        
        MotherWavelet = MotherWaveletNames{kchoice};
    end
set(motherWaveletMenuChoices(1), 'CallBack', @(~,~) selectMotherWaveletMenu(1));
set(motherWaveletMenuChoices(2), 'CallBack', @(~,~) selectMotherWaveletMenu(2));
set(motherWaveletMenuChoices(3), 'CallBack', @(~,~) selectMotherWaveletMenu(3));
set(motherWaveletMenuChoices(4), 'CallBack', @(~,~) selectMotherWaveletMenu(4));
set(motherWaveletMenuChoices(5), 'CallBack', @(~,~) selectMotherWaveletMenu(5));

% frequency scale
FrequencyScaleMenu = uimenu(paramMenu, 'Text', 'Frequency scale');
FrequencyScaleMenuChoices(1) = uimenu(FrequencyScaleMenu, 'Text', 'lin');
FrequencyScaleMenuChoices(2) = uimenu(FrequencyScaleMenu, 'Text', 'log');
FrequencyScaleValues = {'lin', 'log'};
set(FrequencyScaleMenuChoices(find(strcmp(FrequencyScaleValues, FrequencyScale))), 'Checked', 'on');

    function selectFrequencyScaleMenu(kchoice)
        for kchoices = 1:length(FrequencyScaleMenuChoices)
            set(FrequencyScaleMenuChoices(kchoices), 'Checked', 'off');
        end
        set(FrequencyScaleMenuChoices(kchoice), 'Checked', 'on');
        
        FrequencyScale = FrequencyScaleNames{kchoice};
    end
set(FrequencyScaleMenuChoices(1), 'CallBack', @(~,~) selectFrequencyScaleMenu(1));
set(FrequencyScaleMenuChoices(2), 'CallBack', @(~,~) selectFrequencyScaleMenu(2));

% signal channels
signalChannelsMenu = uimenu(paramMenu, 'Text','Set channels', 'Separator', 'on');
    function setSignalChannels()
        signalChannels = getSignalChannels(signalChannels, nbPlotsTot);
        nbPlots = length(signalChannels);
        
        for k_plot = 1:length(WaveletPlot)
            if ismember(k_plot, signalChannels)
                set(WaveletPlot(k_plot), 'LineStyle', '-');
            else
                set(WaveletPlot(k_plot), 'LineStyle', ':');
            end
        end
    end
signalChannelsMenu.Callback = @(~,~) setSignalChannels();

% Xlim
XlimMenu = uimenu(paramMenu, 'Text','Set time limits (Tlim)', 'Separator', 'on');
    function setXlim()
        prompt = {'Enter Tmin :', 'Enter Tmax :'};
        dlgtitle = 'Input time limits (Tlim)';
        dims = [1 35];
        definput = {num2str(XLim(1)), num2str(XLim(2))};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        try
            XLim(1) = str2double(answer{1});
            XLim(2) = str2double(answer{2});
            maxLagCorr = min(maxLagCorr, XLim(2) - XLim(1));
            XLimDisplayObj.updateXLim(XLim);
        catch
        end
    end
set(XlimMenu, 'CallBack', @(~,~) setXlim);

% show XLim
XlimDisplayMenu = uimenu(paramMenu, 'Text','Show Tlim', 'Checked', ShowXLim);
    function setXlimDisplay()
        ShowXLim = ~ShowXLim;
        set(XlimDisplayMenu, 'Checked', ShowXLim);
        XLimDisplayObj.setVisible(ShowXLim, ShowXLim);
    end
set(XlimDisplayMenu, 'CallBack', @(~,~) setXlimDisplay);


% remove mean

removeMeanMenu = uimenu(paramMenu, 'Text','Remove mean', 'Checked', removeMean, 'Separator', 'on');
    function switchRemoveMeanDisplay(status)
        removeMeanMenu.Checked = status;
        removeMean = status;
    end

removeMeanMenu.MenuSelectedFcn = @(~, ~) switchRemoveMeanDisplay(~strcmp(removeMeanMenu.Checked, 'on'));

% ct
CtMenu = uimenu(paramMenu, 'Text', ['Set ct (', num2str(ctEdgeEffects), ')'], 'Separator', 'on');
    function setCt()
        prompt = {'Enter ct :'};
        dlgtitle = 'Input ct';
        dims = [1 35];
        definput = {num2str(ctEdgeEffects)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        try
            ctEdgeEffects = str2num(answer{1});
            set(CtMenu, 'Text', ['Set ct (', num2str(ctEdgeEffects), ')']);
        catch
        end
    end
set(CtMenu, 'CallBack', @(~,~) setCt);

% cf
CfMenu = uimenu(paramMenu, 'Text', ['Set cf (', num2str(cf), ')']);
    function setCf()
        prompt = {'Enter cf :'};
        dlgtitle = 'Input cf';
        dims = [1 35];
        definput = {num2str(cf)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        try
            cf = str2double(answer{1});
            set(CfMenu, 'Text', ['Set cf (', num2str(cf), ')']);
        catch
        end
    end
set(CfMenu, 'CallBack', @(~,~) setCf);

% %zero padding fourier
% zeroPaddingFourierMenu = uimenu(paramMenu, 'Text','Zero padding Fourier');
%     function setZeroPaddingFourier()
%         prompt = {'Enter zero padding Fourier :'};
%         dlgtitle = 'Input zero padding Fourier';
%         dims = [1 35];
%         definput = {num2str(ZeroPaddingFourier)};
%         answer = inputdlg(prompt,dlgtitle,dims,definput);
%         try
%             ZeroPaddingFourier = str2num(answer{1});
%             ZeroPaddingFourier = max(round(ZeroPaddingFourier), 0);
%         catch
%         end
%     end
% set(zeroPaddingFourierMenu, 'CallBack', @(~,~) setZeroPaddingFourier);


%% mode menu

modeMenu = uimenu(fig, 'Text', 'Mode');

%multiSignalMode

multiSignalModeMenu = uimenu(modeMenu, 'Text','Multi signal mode', 'Checked', multiSignalMode);
    function switchMultiSignalModeDisplay(status)
        multiSignalModeMenu.Checked = status;
        multiSignalMode = status;
        set(multiSignalModeOptionsMenu, 'Enable', status);
        if status
            switchMultipleAxesDisplay(false);
            if autocorrelationSVDMode
                switchAutocorrelationModeDisplay(false);
            end
        end
    end

multiSignalModeMenu.MenuSelectedFcn = @(~, ~) switchMultiSignalModeDisplay(~strcmp(multiSignalModeMenu.Checked, 'on'));


% multiSignalMode options

multiSignalModeOptionsMenu = uimenu(modeMenu, 'Text','Multi signal mode options');
multiSignalModeOptionsChoices(1) = uimenu(multiSignalModeOptionsMenu, 'Text', 'sum CWT²');
multiSignalModeOptionsChoices(2) = uimenu(multiSignalModeOptionsMenu, 'Text', 'sum |CWT|²');
    function selectMultiSignalModeOptions(kchoice)
        for kchoices = 1:length(multiSignalModeOptionsChoices)
            set(multiSignalModeOptionsChoices(kchoices), 'Checked', 'off');
        end
        set(multiSignalModeOptionsChoices(kchoice), 'Checked', 'on');
        
        multiSignalModeAbsValue = kchoice == 2;
    end
set(multiSignalModeOptionsChoices(1), 'CallBack', @(~,~) selectMultiSignalModeOptions(1));
set(multiSignalModeOptionsChoices(2), 'CallBack', @(~,~) selectMultiSignalModeOptions(2));

if multiSignalModeAbsValue
    selectMultiSignalModeOptions(2);
else
    selectMultiSignalModeOptions(1);
end


% rdt mode

rdtModeMenu = uimenu(modeMenu, 'Text','Random decrement mode', 'Checked', rdtMode, 'Separator', 'on');
rdtSetMenu = uimenu(modeMenu, 'Text','Set random decrement');

    function switchRdtModeDisplay(status)
        if status
            getXY();
            [Xrdt, Yrdt, axRdt] = RDTmenu(x, y, signalChannels, signalUnit);
            if isempty(Xrdt)
                status = false;
            else
                XLimDisplayObj.addAxesXLimRidges(axRdt, true)
            end
        end
        
        rdtModeMenu.Checked = status;
        rdtMode = status;
        set(rdtSetMenu, 'Enable', status);
        if status
            switchAutocorrelationModeDisplay(false);
        end
    end

    function rdtSetMenuCallback(~,~)
        getXY();
        [Xrdt0, Yrdt0, axRdt] = RDTmenu(x, y, signalChannels, signalUnit);
        if ~isempty(Xrdt0)
            Xrdt = Xrdt0;
            Yrdt = Yrdt0;
            XLimDisplayObj.addAxesXLimRidges(axRdt, true)
        end
    end

rdtModeMenu.MenuSelectedFcn = @(~, ~) switchRdtModeDisplay(~strcmp(rdtModeMenu.Checked, 'on'));
rdtSetMenu.Callback = @rdtSetMenuCallback;

% autocorr mode

autocorrelationModeMenu = uimenu(modeMenu, 'Text','Cross-corr mode', 'Checked', autocorrelationMode, 'Separator', 'on');

autocorrelationParamsMenu = uimenu(modeMenu, 'Text', 'Cross-corr params');

 % plot cross-corr
crosscorrelationDisplayMenu = uimenu(autocorrelationParamsMenu, 'Text', 'Plot cross-corr');

    function displayCrossCorr()
        getXY();
        Dx = (x(1, end) - x(1, 1))/(size(x, 2)-1);
        NmaxLagCorr = floor(maxLagCorr/Dx);
        plotCrossCorr(Dx, y, NmaxLagCorr, 1);
    end
set(crosscorrelationDisplayMenu, 'CallBack', @(~,~) displayCrossCorr);

 % plot autocorr
autocorrelationDisplayMenu = uimenu(autocorrelationParamsMenu, 'Text', 'Plot autocorr');

    function displayAutoCorr()
        getXY();
        Dx = (x(1, end) - x(1, 1))/(size(x, 2)-1);
        NmaxLagCorr = floor(maxLagCorr/Dx);
        plotCrossCorr(Dx, y, NmaxLagCorr, 2);
    end
set(autocorrelationDisplayMenu, 'CallBack', @(~,~) displayAutoCorr);


 % max lag
autocorrelationMaxLagMenu = uimenu(autocorrelationParamsMenu, 'Text', 'Set max lag', 'Separator', 'on');

    function flag = setMaxLagCorr()
        flag = true;
        answer = inputdlg({'Enter max lag:'}, 'Input MaxLag', [1 35], {num2str(maxLagCorr)});
        try
            maxLagCorr = str2double(answer{1});
        catch
            flag = false;
            return
        end
        resetCrossCorr();
        getXY();
        Dx = (x(1, end) - x(1, 1))/(size(x, 2)-1);
        NmaxLagCorr = floor(maxLagCorr/Dx);
        plotCrossCorr(Dx, y, NmaxLagCorr, true);
    end
set(autocorrelationMaxLagMenu, 'CallBack', @(~,~) setMaxLagCorr);

% bias
autocorrelationBiasMenu = uimenu(autocorrelationParamsMenu, 'Text', 'Bias');
autocorrelationBiasValues = {'biased', 'unbiased'};
autocorrelationBiasMenuChoices(1) = uimenu(autocorrelationBiasMenu, 'Text', 'Biased');
autocorrelationBiasMenuChoices(2) = uimenu(autocorrelationBiasMenu, 'Text', 'Unbiased');
set(autocorrelationBiasMenuChoices(find(strcmp(autocorrelationBiasValues, autocorrelationBias))), 'Checked', 'on');

    function selectAutocorrelationBiasMenu(kchoice)
        for kchoices = 1:length(autocorrelationBiasMenuChoices)
            set(autocorrelationBiasMenuChoices(kchoices), 'Checked', 'off');
        end
        set(autocorrelationBiasMenuChoices(kchoice), 'Checked', 'on');
        
        autocorrelationBias = autocorrelationBiasValues{kchoice};
        resetCrossCorr();
        getXY();
        Dx = (x(1, end) - x(1, 1))/(size(x, 2)-1);
        NmaxLagCorr = floor(maxLagCorr/Dx);
        plotCrossCorr(Dx, y, NmaxLagCorr, true);
    end
set(autocorrelationBiasMenuChoices(1), 'CallBack', @(~,~) selectAutocorrelationBiasMenu(1));
set(autocorrelationBiasMenuChoices(2), 'CallBack', @(~,~) selectAutocorrelationBiasMenu(2));


 % set nb of SV

AutocorrelationNsvdMenu = uimenu(autocorrelationParamsMenu, 'Text', sprintf('Set nb of SV (%u)', autocorrelationNsvd));

    function AutocorrelationNsvdMenuCallback(~,~)
        while true
            answer = inputdlg({'Enter nb. of singular values:'}, 'Input nb of SV', [1 35], {num2str(autocorrelationNsvd)});
            try
                autocorrelationNsvd0 = str2double(answer{1});
                if isnan(autocorrelationNsvd0) || mod(autocorrelationNsvd0, 1) ~= 0 ||...
                        autocorrelationNsvd0 <= 0 || autocorrelationNsvd0 > nbPlots
                    figError = errordlg('Incorrect input', 'Error', 'modal');
                    waitfor(figError);
                else
                    break
                end
            catch
                return
            end
        end
        autocorrelationNsvd = autocorrelationNsvd0;
        set(AutocorrelationNsvdMenu, 'Text', sprintf('Set nb of SV (%u)', autocorrelationNsvd));
    end
AutocorrelationNsvdMenu.MenuSelectedFcn = @AutocorrelationNsvdMenuCallback;

 % svd modes

autocorrelationSVDMenu = uimenu(autocorrelationParamsMenu, 'Text', 'SVD mode CWT', 'Separato', 'on');
if autocorrelationSVDMode
    set('autocorrelationSVDMenu', 'Checked', 'on');
end
    function autocorrSVDCallback(~,~)
        if autocorrelationSVDMode
            set(autocorrelationSVDMenu, 'Checked', 'off');
            autocorrelationSVDMode = false;
        else
            set(autocorrelationSVDMenu, 'Checked', 'on');
            autocorrelationSVDMode = true;
        end
    end
autocorrelationSVDMenu.Callback = @autocorrSVDCallback;


autocorrelationFourierMenu = uimenu(autocorrelationParamsMenu, 'Text', 'SVD mode Fourier', 'Checked', autocorrelationFourierSVDMode);

    function switchAutocorrelationModeDisplay(status)
        if status
            if setMaxLagCorr()
                switchRdtModeDisplay(false);
                switchMultipleAxesDisplay(false);
                if autocorrelationSVDMode
                    switchMultiSignalModeDisplay(false);
                end
            else
                status = false;
            end
        end
        autocorrelationModeMenu.Checked = status;
        autocorrelationMode = status;
        set(autocorrelationParamsMenu, 'Enable', status);
        switchAutocorrelationFourierMenu(autocorrelationFourierSVDMode);
    end
autocorrelationModeMenu.MenuSelectedFcn = @(~, ~) switchAutocorrelationModeDisplay(~strcmp(autocorrelationModeMenu.Checked, 'on'));

    function switchAutocorrelationFourierMenu(flag)
        if nargin == 0
            flag = ~autocorrelationFourierSVDMode;
        end
        
        autocorrelationFourierSVDMode = flag;
            set(autocorrelationFourierMenu, 'Checked', flag);
        if flag && autocorrelationMode
            set(FourierAveragingMenu, 'Enable', 'off');
        else
            set(FourierAveragingMenu, 'Enable', 'on');
        end
        
    end
set(autocorrelationFourierMenu, 'CallBack', @(~,~) switchAutocorrelationFourierMenu());

%set

switchMultiSignalModeDisplay(multiSignalMode);
switchRdtModeDisplay(rdtMode);
% switchAutocorrelationModeDisplay(autocorrelationMode); % after FourierAveragingMenu


%% ridge menu

ridgeMenu = uimenu(fig,'Text','Ridges');

%autres valeurs par défault
freqRidgeName = 'freq';
phaseRidgeName = 'pha2';
dampingRidgeName = 'damping3';

freqRidgeNames = {'freq', 'freq2'};
phaseRidgeNames = {'pha', 'pha2'};
dampingRidgeNames = {'damping', 'damping2', 'damping3'};
FrequencyScaleNames = {'lin', 'log'};

%freq
freqMenu = uimenu(ridgeMenu, 'Text','Frequency');
freqMenuChoices(1) = uimenu(freqMenu, 'Text', 'max module', 'Checked' ,'on');
freqMenuChoices(2) = uimenu(freqMenu, 'Text', 'phase derivative');
    function selectFreqMenu(kchoice)
        for kchoices = 1:length(freqMenuChoices)
            set(freqMenuChoices(kchoices), 'Checked', 'off');
        end
        set(freqMenuChoices(kchoice), 'Checked', 'on');
        
        freqRidgeName = freqRidgeNames{kchoice};
    end
set(freqMenuChoices(1), 'CallBack', @(~,~) selectFreqMenu(1));
set(freqMenuChoices(2), 'CallBack', @(~,~) selectFreqMenu(2));

%phase
phaseMenu = uimenu(ridgeMenu, 'Text','Phase');
phaseMenuChoices(1) = uimenu(phaseMenu, 'Text', 'bounded');
phaseMenuChoices(2) = uimenu(phaseMenu, 'Text', 'continuous', 'Checked' ,'on');
    function selectPhaseMenu(kchoice)
        for kchoices = 1:length(phaseMenuChoices)
            set(phaseMenuChoices(kchoices), 'Checked', 'off');
        end
        set(phaseMenuChoices(kchoice), 'Checked', 'on');
        
        phaseRidgeName = phaseRidgeNames{kchoice};
    end
set(phaseMenuChoices(1), 'CallBack', @(~,~) selectPhaseMenu(1));
set(phaseMenuChoices(2), 'CallBack', @(~,~) selectPhaseMenu(2));

% damping
dampingMenu = uimenu(ridgeMenu, 'Text','Damping');
dampingMenuChoices(1) = uimenu(dampingMenu, 'Text', 'lambda');
dampingMenuChoices(2) = uimenu(dampingMenu, 'Text', 'lambda/omega_d');
dampingMenuChoices(3) = uimenu(dampingMenu, 'Text', 'lambda/omega_n', 'Checked' ,'on');
    function selectDampingMenu(kchoice)
        for kchoices = 1:length(dampingMenuChoices)
            set(dampingMenuChoices(kchoices), 'Checked', 'off');
        end
        set(dampingMenuChoices(kchoice), 'Checked', 'on');
        
        dampingRidgeName = dampingRidgeNames{kchoice};
    end
set(dampingMenuChoices(1), 'CallBack', @(~,~) selectDampingMenu(1));
set(dampingMenuChoices(2), 'CallBack', @(~,~) selectDampingMenu(2));
set(dampingMenuChoices(3), 'CallBack', @(~,~) selectDampingMenu(3));

% Xlim
XlimRidgeMenu = uimenu(ridgeMenu, 'Text','Set time limits (Tlim) ridge', 'Separator', 'on');
    function setXlimRidge()
        prompt = {'Enter Tmin ridge :', 'Enter Tmax ridge :'};
        dlgtitle = 'Input time limits (Tlim) ridge';
        dims = [1 35];
        definput = {num2str(XLimRidge(1)), num2str(XLimRidge(2))};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        try
            XLimRidge(1) = str2double(answer{1});
            XLimRidge(2) = str2double(answer{2});
            XLimDisplayObj.updateXLimRidges(XLimRidge);
        catch
        end
    end
set(XlimRidgeMenu, 'CallBack', @(~,~) setXlimRidge);

% ct
CtRidgeMenu = uimenu(ridgeMenu, 'Text', ['Set ct ridge (', num2str(ctRidge), ')']);
    function setCtRidge()
        prompt = {'Enter ct ridge:'};
        dlgtitle = 'Input ct ridge';
        dims = [1 35];
        definput = {num2str(ctRidge)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        try
            ctRidge = str2num(answer{1});
            set( CtRidgeMenu, 'Text', ['Set ct ridge (', num2str(ctRidge), ')']);
        catch
        end
    end
set(CtRidgeMenu, 'CallBack', @(~,~) setCtRidge);

% stop when incrreasing
StopRidgeWhenIncreasingMenu = uimenu(ridgeMenu, 'Text', 'Stop ridge when increasing',...
    'Separator', 'on', 'Checked', StopRidgeWhenIncreasing);
    function toggleStopRidgeWhenIncreasing()
        StopRidgeWhenIncreasing = ~StopRidgeWhenIncreasing;
        set(StopRidgeWhenIncreasingMenu, 'Checked', StopRidgeWhenIncreasing);
    end
StopRidgeWhenIncreasingMenu.Callback = @(~,~) toggleStopRidgeWhenIncreasing();

% threshold
ThresholdRidgeMenu = uimenu(ridgeMenu, 'Text','Set ridge threshold', 'Separator', 'on');
    function setThresholdRidge()
        prompt = {['Enter ridge threshold [', signalUnit, ']:']};
        dlgtitle = 'Input ridge threshold';
        dims = [1 35];
        definput = {sprintf('%.10e', RidgeMinModu)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        try
            RidgeMinModu = eval(answer{1});
        catch
        end
    end
set(ThresholdRidgeMenu, 'CallBack', @(~,~) setThresholdRidge);

% compute average noise
AverageNoiseMenu = uimenu(ridgeMenu, 'Text','Get average |CWT|');
    function getAverageNoise()
        % input
        QAverage = eval(get(editQ, 'String'));
        
        % menu
        [fAverage, ignoreXLim, typeOfAverage, figAverageMenu] = getAverageMenu(WvltScale);
        
        if isempty(fAverage) || isnan(fAverage)
            try
                delete(figAverageMenu);
            catch
            end
            return
        end
        
        % compute average
        if ignoreXLim
            XminOld = XLim(1); XmaxOld = XLim(2);
            XLim(1) = -inf; XLim(2) = inf;
            getXY();
            XLim(1) = XminOld; XLim(2) = XmaxOld;
        else
            getXY();
        end
        
        if multiSignalMode
            wavelet = 0;
            for kPlot = 1:nbPlots
                if ~multiSignalModeAbsValue
                    wavelet = wavelet + WvltComp(x(kPlot,:), y(kPlot,:), fAverage, QAverage,...
                        'MotherWavelet', MotherWavelet).^2;
                else
                    wavelet = wavelet + abs(WvltComp(x(kPlot,:), y(kPlot,:), fAverage, QAverage,...
                        'MotherWavelet', MotherWavelet)).^2;
                end
            end
        else
            wavelet = [];
            for kPlot = 1:nbPlots
                wavelet = [wavelet; WvltComp(x(kPlot,:), y(kPlot,:), fAverage, QAverage,...
                    'MotherWavelet', MotherWavelet)];
            end
        end
        
        % edge effects
        [~, DeltaT] = FTpsi_DeltaT(QAverage, MotherWavelet);
        deltaT = DeltaT(fAverage);
        deltaT = ctEdgeEffects .* [deltaT, deltaT];
        x1 = x(1, :);
        wavelet = wavelet(:, (x1 >= x1(1) + deltaT(1)) & (x1 <= x1(end) - deltaT(2)));
        
        % compute
        switch typeOfAverage
            case 'lin'
                AverageNoise = mean(abs(wavelet), 2);
            case 'log'
                AverageNoise = exp(mean(log(abs(wavelet)), 2));
        end
        
        if multiSignalMode
            AverageNoise = sqrt(AverageNoise);
        end
        
        % display results
        try
            delete(figAverageMenu);
        catch
        end
        
        dispAverageMenu(AverageNoise, fAverage, multiSignalMode, typeOfAverage, ignoreXLim, signalUnit);
        
    end
set(AverageNoiseMenu, 'CallBack', @(~,~) getAverageNoise);


%multipleAxesDisplay

multipleAxesDisplayMenu = uimenu(ridgeMenu, 'Text','Multiple axes',...
    'Checked', multipleAxesDisplay, 'Separator', 'on');
    function switchMultipleAxesDisplay(status)
        multipleAxesDisplayMenu.Checked = status;
        setMultipleAxesDisplay(status);
        if status
            switchMultiSignalModeDisplay(false);
            switchAutocorrelationModeDisplay(false)
        end
    end

multipleAxesDisplayMenu.MenuSelectedFcn = @(~, ~) switchMultipleAxesDisplay(~strcmp(multipleAxesDisplayMenu.Checked, 'on'));



%% fourier menu

FourierScaleNames = {'lin', 'squared', 'log', 'spectral density (lin)', 'spectral density (log)', 'phase'};
FourierWindowNames = {'rectangular', 'hamming', 'exponential'};


fourierMenu = uimenu(fig,'Text','Fourier');

% fourier scale
FourierScaleMenu = uimenu(fourierMenu, 'Text','Fourier Scale');
FourierScaleMenuChoices(1) = uimenu(FourierScaleMenu, 'Text', 'lin');
FourierScaleMenuChoices(2) = uimenu(FourierScaleMenu, 'Text', 'squared');
FourierScaleMenuChoices(3) = uimenu(FourierScaleMenu, 'Text', 'log');
FourierScaleMenuChoices(4) = uimenu(FourierScaleMenu, 'Text', 'spectral density (lin)');
FourierScaleMenuChoices(5) = uimenu(FourierScaleMenu, 'Text', 'spectral density (log)');
FourierScaleMenuChoices(6) = uimenu(FourierScaleMenu, 'Text', 'phase');
set(FourierScaleMenuChoices(find(strcmp(FourierScaleNames, FourierScale))), 'Checked', 'on');

    function selectFourierScaleMenu(kchoice)
        for kchoices = 1:length(FourierScaleMenuChoices)
            set(FourierScaleMenuChoices(kchoices), 'Checked', 'off');
        end
        set(FourierScaleMenuChoices(kchoice), 'Checked', 'on');
        
        FourierScale = FourierScaleNames{kchoice};
    end
set(FourierScaleMenuChoices(1), 'CallBack', @(~,~) selectFourierScaleMenu(1));
set(FourierScaleMenuChoices(2), 'CallBack', @(~,~) selectFourierScaleMenu(2));
set(FourierScaleMenuChoices(3), 'CallBack', @(~,~) selectFourierScaleMenu(3));
set(FourierScaleMenuChoices(4), 'CallBack', @(~,~) selectFourierScaleMenu(4));
set(FourierScaleMenuChoices(5), 'CallBack', @(~,~) selectFourierScaleMenu(5));
set(FourierScaleMenuChoices(6), 'CallBack', @(~,~) selectFourierScaleMenu(6));


% fourier window
FourierWindowMenu = uimenu(fourierMenu, 'Text','Window', 'Separator', 'on');
FourierWindowMenuChoices(1) = uimenu(FourierWindowMenu, 'Text', 'rectangular');
FourierWindowMenuChoices(2) = uimenu(FourierWindowMenu, 'Text', 'hamming');
FourierWindowMenuChoices(3) = uimenu(FourierWindowMenu, 'Text', 'exponential');
set(FourierWindowMenuChoices(find(strcmp(FourierWindowNames, FourierWindow))), 'Checked', 'on');

    function selectFourierWindowMenu(kchoice)
        for kchoices = 1:length(FourierWindowMenuChoices)
            set(FourierWindowMenuChoices(kchoices), 'Checked', 'off');
        end
        set(FourierWindowMenuChoices(kchoice), 'Checked', 'on');
        
        FourierWindow = FourierWindowNames{kchoice};
        
        switch FourierWindow
            case 'exponential'
                while true
                    answer = inputdlg({'Enter \tau [s]:'}, 'Input exponential window params', [1 35], {''});
                    try
                        FourierWindowParams = str2double(answer{1});
                        break
                    catch
                    end
                end
        end
    end
set(FourierWindowMenuChoices(1), 'CallBack', @(~,~) selectFourierWindowMenu(1));
set(FourierWindowMenuChoices(2), 'CallBack', @(~,~) selectFourierWindowMenu(2));
set(FourierWindowMenuChoices(3), 'CallBack', @(~,~) selectFourierWindowMenu(3));

% Fourier Averaging

FourierAveragingMenu = uimenu(fourierMenu, 'Text', 'Averaging', 'Checked', FourierAveraging, 'Separator', 'on');
FourierAveragingNbMenu = uimenu(fourierMenu, 'Text', sprintf('Set averaging (%u)', FourierAveragingNb));

    function updateFourierAveragingMenu(~,~)
        FourierAveraging = ~FourierAveraging;
        set(FourierAveragingMenu, 'Checked', FourierAveraging);
        if FourierAveraging
            set(FourierAveragingNbMenu, 'Enable', 'on');
            FourierAveragingNbMenuCallback;
        else
            set(FourierAveragingNbMenu, 'Enable', 'off');
        end
    end

    function FourierAveragingNbMenuCallback(~,~)
        dlgtitle = 'Input averaging number';
        prompt = {'Divide signal into N subsignals:'};
        dims = [1 35];
        definput = {num2str(FourierAveragingNb)};
        while true
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            try
                FourierAveragingNb0 = str2double(answer{1});
            catch
                return
            end
            
            if mod(FourierAveragingNb0, 1) ~= 0 || FourierAveragingNb0 < 1
                errorFig = errordlg('Incorrect input', 'Error', 'modal');
                waitfor(errorFig);
            else
                break
            end
        end
        
        FourierAveragingNb = FourierAveragingNb0;
        set(FourierAveragingNbMenu, 'Text', sprintf('Set averaging (%u)', FourierAveragingNb));
    end

FourierAveragingMenu.MenuSelectedFcn = @updateFourierAveragingMenu;
FourierAveragingNbMenu.MenuSelectedFcn = @FourierAveragingNbMenuCallback;

FourierAveraging = ~FourierAveraging;
updateFourierAveragingMenu(FourierAveraging);

switchAutocorrelationModeDisplay(autocorrelationMode);


%% shock menu

shocksMenu = uimenu(fig, 'Text', 'Shocks');

% shocks detection
getShocksMenu = uimenu(shocksMenu, 'Text', 'Shocks detection menu');

    function getShocksMenuCallback()
        % evaluation des parametres
        fminShock = eval(get(editfmin, 'String'));
        fmaxShock = eval(get(editfmax, 'String'));
        NbFreqShock = eval(get(editNbFreq, 'String'));
        QShock = eval(get(editQ, 'String'));
        getXY();
        
        % shocks menu
        shockDetectionMenu(x, y, signalChannels,...
            QShock, MotherWavelet, ctEdgeEffects, QShock, MotherWavelet,...
            [fminShock, fmaxShock], NbFreqShock, FrequencyScale, WvltScale,...
            [fminShock, fmaxShock], NbFreqShock, FrequencyScale, WvltScale,...
            multiSignalMode, multiSignalMode);
    end

getShocksMenu.MenuSelectedFcn = @(~, ~) getShocksMenuCallback();


%% tools menu

toolsMenu = uimenu(fig, 'Text', 'Tools');

% bounds Q

QboundsMenu = uimenu(toolsMenu,'Text','Bounds Q');

QboundsMenu.MenuSelectedFcn = @QboundsMenuCallback;
    function QboundsMenuCallback(~,~)
        BoundsQMenu(XLim, XLimRidge, ctEdgeEffects, ctRidge, cf, MotherWavelet);
    end

% filtering

filteringMenu = uimenu(toolsMenu,'Text','Filtering');

filteringMenu.MenuSelectedFcn = @(~, ~) LinearFilterMain(WaveletPlot);

% regressions

regMenu = uimenu(toolsMenu,'Text','Regression');

regConstant = uimenu(regMenu, 'Text','Constant');
regExp = uimenu(regMenu, 'Text','Exp');

regConstant.MenuSelectedFcn = @(~, ~) RegressionMenu('Equation', 'c', 'Param', 'c', 'Param0', 0);
regExp.MenuSelectedFcn = @(~, ~) RegressionMenu('Equation', 'a*exp(-lambda*x)', 'Param', 'a  lambda',...
    'Param0', [1 1], 'Fit', 'log(y)');


% averaging
averagingMenu = uimenu(toolsMenu, 'Text', 'Averaging');

averagingMenu.MenuSelectedFcn = @(~, ~) ScatterAveragingMenu();


% plot extract

plotExtractMenu = uimenu(toolsMenu,'Text','Plot Extract');

    function plotExtractCallback(~, ~)
        prompt = {'Enter x-axis var name:', 'Enter y-axis var name:'};
        dlgtitle = 'Input var names';
        dims = [1 35];
        definput = {'', ''};
        answer = inputdlg(prompt, dlgtitle, dims, definput);
        
        if isempty(answer) % cancel button
            return
        end
        
        try
            evalin('base', ['[', answer{1}, ',', answer{2}, ']=PlotExtract();']);
        catch e
            warning(e.message);
        end
    end

plotExtractMenu.MenuSelectedFcn = @plotExtractCallback;


% audio
audioMenu = uimenu(toolsMenu, 'Text', 'Audio');

    function audioMenuCallback()
        audioMain(WaveletPlot)
    end

audioMenu.MenuSelectedFcn = @(~, ~) audioMenuCallback();


%% help menu

helpMenu = uimenu(fig, 'Text', 'Help');

helpDocMenu = uimenu(helpMenu, 'Text', 'Documentation');
helpDocMenu.Callback = @(~, ~) documentation_waveletmenu(false);

helpUpdateMenu = uimenu(helpMenu, 'Text', 'Update', 'Separator', 'on');
helpUpdateMenu.Callback = @(~, ~) updateCWT();

if CheckForUpdates && ~checkUpdate()
    helpMenu.ForegroundColor = 'red';
    helpMenu.Text = [helpMenu.Text, '*'];
    helpUpdateMenu.ForegroundColor = 'red';
    helpUpdateMenu.Text = [helpUpdateMenu.Text, '*'];
end

%% info

% parameters panel
strfmin.Tooltip = 'minimum frequency [Hz] (CWT & Fourier)';
strfmax.Tooltip = 'maximum frequency [Hz] (CWT & Fourier)';
strNbFreq.Tooltip = 'number of discrete frequencies (CWT)';
strQ.Tooltip = ['quality factor', newline, 'see: Le & Argoul 2004, Erlicher & Argoul 2007'];
strmaxR.Tooltip = 'maximum number of ridges';
strPR.Tooltip = 'maximum number of parallel (equal time) ridges';
strSL.Tooltip = ['ridge stopping condition', newline, 'df/dt < [max ridge slope]*f²'];

% plots panel
checkboxModule.Tooltip = '|CWT(t, f)|';
checkboxPhase.Tooltip = 'arg CWT(t, f)';
checkboxTimeAmplPlot.Tooltip = 'ridge plot : time X amplitude (directly on original axes)';
checkboxTimeAmpl.Tooltip = 'ridge plot : time X amplitude';
checkboxTimeFreq.Tooltip = 'ridge plot : time X frequency';
checkboxTimeDamp.Tooltip = 'ridge plot : time X damping';
checkboxAmplFreq.Tooltip = 'ridge plot : amplitude X frequency';
checkboxAmplDamp.Tooltip = 'ridge plot : amplitude X damping';
checkboxTimePhase.Tooltip = 'ridge plot : time X phase';

% mode shapes panel
strShapesInstantaneous.Tooltip = 'instantaneous mode shapes';
strxAxisShapes.Tooltip = 'instantaneous mode shapes';
stryAxisShapes.Tooltip = ['instantaneous mode shapes', newline, 'shaping: phi*phi^T = 1 (phi in C^n)'];
checkboxRealShapes.Tooltip = ['instantaneous mode shapes: real parts', newline, 'shaping: phi*phi^T = 1 (phi in C^n)'];
checkboxImagShapes.Tooltip = ['instantaneous mode shapes: imaginary parts', newline, 'shaping: phi*phi^T = 1 (phi in C^n)'];
checkboxModuleShapes.Tooltip = ['instantaneous mode shapes: modules', newline, 'shaping: phi*phi^T = 1 (phi in C^n)'];
checkboxPhaseShapes.Tooltip = ['instantaneous mode shapes: angles', newline, 'shaping: phi*phi^T = 1 (phi in C^n)'];
checkboxAmplShapes.Tooltip = ['instantaneous mode shapes: amplitude', newline, 'shaping: phi*phi^T = 1 (phi in C^n)'];
strShapesMean.Tooltip = ['averaged mode shapes', newline, 'shaping before average: phi*phi^T = 1 (phi in C^n)'];
weightOptionShapesMean.Tooltip = 'average weighted by ridge amplitude';
checkboxRealShapesMean.Tooltip = ['averaged mode shapes (real part) plot', newline, 'shaping before average: phi*phi^T = 1 (phi in C^n)',...
    newline, 'customizable: ''RealShapePlot'' option'];
checkboxComplexShapesMean.Tooltip = ['averaged mode shapes (real and imaginary parts) plot',...
    newline, 'shaping before average: phi*phi^T = 1 (phi in C^n)', newline, 'customizable: ''ComplexShapePlot'' option'];
checkboxDispShapesMean.Tooltip = ['averaged mode shapes (real and imaginary parts) in Command Window',...
    newline, 'shaping before average: phi*phi^T = 1 (phi in C^n)'];
checkboxDispFreqsMean.Tooltip = 'averaged frequencies in Command Window';
checkboxAmplRegMean.Tooltip = 'linear regression on amplitude log';



%%

    function getXY()
        X = getX();
        Y = getY();
        
        % test X Y size
        if ~isequal(size(X), size(Y))
            error('input matrices size error');
        end
        
        % test X identical lines
        for line = 2:size(X, 1)
            if ~isequal(X(line, :), X(1, :))
                error('time arrays must be the same for every plot');
            end
        end
        
        % test time step
        if any( abs(diff(X(1, :)) / mean(diff(X(1, :))) - 1) > 1e-4)
            error('non-constant time step');
        end
        
        
        for line = 1:size(X, 1)
            Y2(line, :) = Y(line, X(line, :)>=XLim(1) & X(line, :)<=XLim(2));
            X2(line, :) = X(line, X(line, :)>=XLim(1) & X(line, :)<=XLim(2));
        end
        
        if removeMean
            Y2 = Y2 - mean(Y2, 2);
        end
        
        X2 = X2(signalChannels, :);
        Y2 = Y2(signalChannels, :);
        
        if ~isequal(x, X2) || ~isequal(y, Y2)
            resetCrossCorr();
        end
        
        x = X2;
        y = Y2;
    end

    function show()
        % changement du curseur
        set(fig, 'pointer', 'watch');
        drawnow;
        
        % evaluation des parametres
        fminNew = eval(get(editfmin, 'String'));
        fmaxNew = eval(get(editfmax, 'String'));
        NbFreqNew = eval(get(editNbFreq, 'String'));
        QNew = eval(get(editQ, 'String'));
        if fminNew ~= fmin || fmaxNew ~= fmax || NbFreqNew ~= NbFreq || QNew ~= Q
            resetCrossCorr();
        end
        fmin = fminNew; fmax = fmaxNew; NbFreq = NbFreqNew; Q = QNew;
        
        switch FrequencyScale
            case 'lin'
                WvltFreqs = linspace(fmin, fmax, NbFreq);
            case 'log'
                WvltFreqs = logspace(log10(fmin), log10(fmax), NbFreq);
        end
        
        maxR = eval(get(editmaxR, 'String')); %nombre max de ridges
        PR = eval(get(editPR, 'String')); %nombre max de ridges parallèles
        slopeRidge = eval(get(editSL, 'String')); %max slope ridge
        getXY();
        Dx = (x(1, end) - x(1, 1))/(size(x, 2)-1);
        
        % cross corr
        if autocorrelationMode
            NmaxLagCorr = floor(maxLagCorr/Dx);
            plotCrossCorr(Dx, y, NmaxLagCorr);
            if autocorrelationSVDMode && isempty(SVry)
                [SVry, SVvectry] = svdCWT(tRy, Ry, WvltFreqs, Q, autocorrelationNsvd);
            end
        end
        
        %% plot de la transformee
        if checkboxGeneral.Value
            if ~multiSignalMode
                for kPlot = 1:nbPlots
                    WvltPlot(x(kPlot,:), y(kPlot,:), WvltFreqs, Q, 'ctEdgeEffects', ctEdgeEffects, 'FreqScale', FrequencyScale);
                end
            end
        end
        
        if checkboxModule.Value || checkboxPhase.Value
            if (~multiSignalMode && ~autocorrelationMode) ||...
                (~multiSignalMode && autocorrelationMode && ~autocorrelationSVDMode)
                for kPlot = 1:nbPlots
                    if ~autocorrelationMode && ~rdtMode
                        xCWTplot = x(kPlot,:);
                        yCWTplot = y(kPlot,:);
                        titleCWTPlot = ['channel ', num2str(signalChannels(kPlot)), ' ; Q=', num2str(Q),' ; scale:', WvltScale];
                    elseif autocorrelationMode && ~autocorrelationSVDMode
                        xCWTplot = tRy;
                        yCWTplot = Ry(kPlot,kPlot,:);
                        yCWTplot = transpose(yCWTplot(:));
                        titleCWTPlot = ['R', num2str(signalChannels(kPlot)), num2str(signalChannels(kPlot)), ' ; Q=', num2str(Q),' ; scale:', WvltScale];
                    elseif rdtMode
                        xCWTplot = Xrdt;
                        yCWTplot = Yrdt(kPlot,:);
                        titleCWTPlot = ['RDT channel ', num2str(signalChannels(kPlot)), num2str(signalChannels(kPlot)), ' ; Q=', num2str(Q),' ; scale:', WvltScale];
                    end
                    wavelet = WvltComp(xCWTplot, yCWTplot, WvltFreqs, Q, 'MotherWavelet', MotherWavelet);
                    
                    if checkboxModule.Value
                        WvltPlot2(xCWTplot, WvltFreqs, wavelet, 'module', Q, ctEdgeEffects, MotherWavelet,...
                            WvltScale, titleCWTPlot, wvltAxesTitle, FrequencyScale);
                    end
                    if checkboxPhase.Value
                        WvltPlot2(xCWTplot, WvltFreqs, wavelet, 'phase', Q, ctEdgeEffects, MotherWavelet,...
                            WvltScale, titleCWTPlot, wvltAxesTitle, FrequencyScale);
                    end
                end
                
            elseif multiSignalMode
                wavelet = 0; % calcul de la somme des carrés de transformées
                for kPlot = 1:nbPlots
                    if ~autocorrelationMode && ~rdtMode
                        xCWTplot = x(kPlot,:);
                        yCWTplot = y(kPlot,:);
                    elseif autocorrelationMode && ~autocorrelationSVDMode
                        xCWTplot = tRy;
                        yCWTplot = Ry(kPlot,kPlot,:);
                        yCWTplot = transpose(yCWTplot(:));
                    elseif rdtMode
                        xCWTplot = Xrdt;
                        yCWTplot = Yrdt(kPlot,:);
                    end
                    
                    if ~multiSignalModeAbsValue
                        wavelet = wavelet +...
                            WvltComp(xCWTplot, yCWTplot, WvltFreqs, Q, 'MotherWavelet', MotherWavelet).^2;
                    else
                        wavelet = wavelet +...
                            abs(WvltComp(xCWTplot, yCWTplot, WvltFreqs, Q, 'MotherWavelet', MotherWavelet)).^2;
                    end
                end
                
                if ~autocorrelationMode && ~rdtMode
                    titleCWTPlot = ['sum_wvlt^2 (all selected channels) ; Q=', num2str(Q),' ; scale: ', WvltScale];
                elseif autocorrelationMode && ~autocorrelationSVDMode
                    titleCWTPlot = ['sum_wvlt^2 (all Rxx) ; Q=', num2str(Q),' ; scale: ', WvltScale];
                elseif rdtMode
                    titleCWTPlot = ['sum_wvlt^2 (all RDT channels) ; Q=', num2str(Q),' ; scale: ', WvltScale];
                end
                if checkboxModule.Value
                    WvltPlot2(xCWTplot, WvltFreqs, wavelet,...
                        'module', Q, ctEdgeEffects, MotherWavelet, WvltScale,...
                        titleCWTPlot, wvltAxesTitle, FrequencyScale);
                end
                if checkboxPhase.Value
                    WvltPlot2(xCWTplot, WvltFreqs, wavelet,...
                        'phase', Q, ctEdgeEffects, MotherWavelet, WvltScale,...
                        titleCWTPlot, wvltAxesTitle, FrequencyScale);
                end
                
            elseif autocorrelationMode && autocorrelationSVDMode
                for ksvd = 1:autocorrelationNsvd
                    titleCWTPlot = ['xcorr->CWT->SVD (all selected channels) ; sing. value ', num2str(ksvd),...
                        ' ; Q=', num2str(Q),' ; scale: ', WvltScale];
                    if checkboxModule.Value
                        WvltPlot2(tRy, WvltFreqs, SVry{ksvd},...
                            'module', Q, ctEdgeEffects, MotherWavelet, WvltScale,...
                            titleCWTPlot, wvltAxesTitle, FrequencyScale);
                    end
                    if checkboxPhase.Value
                        WvltPlot2(tRy, WvltFreqs, SVry{ksvd},...
                            'phase', Q, ctEdgeEffects, MotherWavelet, WvltScale,...
                            titleCWTPlot, wvltAxesTitle, FrequencyScale);
                    end
                end
            end
        end
        
        %% calcul des ridges
        if any([Checkboxs2.Value]) || (exist('checkboxTimeAmplPlot', 'var') && checkboxTimeAmplPlot.Value)
            
            if (~multiSignalMode && ~autocorrelationMode) ||...
                    (~multiSignalMode && autocorrelationMode && ~autocorrelationSVDMode)
                ridges = cell(1, nbPlots);
                for kPlot = 1:nbPlots
                    if ~autocorrelationMode && ~rdtMode
                        xRidgePlot = x(kPlot,:);
                        yRidgePlot = y(kPlot,:);
                    elseif autocorrelationMode && ~autocorrelationSVDMode
                        xRidgePlot = tRy;
                        yRidgePlot = Ry(kPlot,kPlot,:);
                        yRidgePlot = transpose(yRidgePlot(:));
                    elseif rdtMode
                        xRidgePlot = Xrdt;
                        yRidgePlot = Yrdt(kPlot,:);
                    end
                    ridges{kPlot} = RidgeExtract(xRidgePlot, yRidgePlot, Q, fmin, fmax, NbFreq,...
                        'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu,...
                        'StopWhenIncreasing', StopRidgeWhenIncreasing,...
                        'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects, 'MaxSlopeRidge', slopeRidge,...
                        'MotherWavelet', MotherWavelet, 'XLimRidge', XLimRidge, 'ctRidge', ctRidge,...
                        'FrequencyScale', FrequencyScale, 'SquaredCWT', multiSignalMode);
                end
                
            elseif multiSignalMode
                wavelet = 0;
                for kPlot = 1:nbPlots
                    if ~autocorrelationMode && ~rdtMode
                        xRidgePlot = x(kPlot,:);
                        yRidgePlot = y(kPlot,:);
                    elseif autocorrelationMode && ~autocorrelationSVDMode
                        xRidgePlot = tRy;
                        yRidgePlot = Ry(kPlot,kPlot,:);
                        yRidgePlot = transpose(yRidgePlot(:));
                    elseif rdtMode
                        xRidgePlot = Xrdt;
                        yRidgePlot = Yrdt(kPlot,:);
                    end
                    
                    if ~multiSignalModeAbsValue
                        wavelet = wavelet +...
                            WvltComp(xRidgePlot, yRidgePlot, WvltFreqs, Q, 'MotherWavelet', MotherWavelet).^2;
                    else
                        wavelet = wavelet +...
                            abs(WvltComp(xRidgePlot, yRidgePlot, WvltFreqs, Q, 'MotherWavelet', MotherWavelet)).^2;
                    end
                end
                
                ridges = cell(1, 1);
                ridges{1} = RidgeExtract(xRidgePlot, nan, Q, fmin, fmax, NbFreq,...
                    'Wavelet', wavelet,...
                    'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu,...
                    'StopWhenIncreasing', StopRidgeWhenIncreasing,...
                    'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects, 'MaxSlopeRidge', slopeRidge,...
                    'MotherWavelet', MotherWavelet, 'XLimRidge', XLimRidge, 'ctRidge', ctRidge,...
                    'FrequencyScale', FrequencyScale, 'SquaredCWT', multiSignalMode);
            
            elseif autocorrelationMode && autocorrelationSVDMode
                ridges = cell(1, autocorrelationNsvd);
                for ksvd = 1:autocorrelationNsvd
                    ridges{ksvd} = RidgeExtract(tRy, nan, Q, fmin, fmax, NbFreq, 'Wavelet', SVry{ksvd},...
                        'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu,...
                    'StopWhenIncreasing', StopRidgeWhenIncreasing,...
                        'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects, 'MaxSlopeRidge', slopeRidge,...
                        'MotherWavelet', MotherWavelet, 'XLimRidge', XLimRidge, 'ctRidge', ctRidge,...
                        'FrequencyScale', FrequencyScale, 'SquaredCWT', multiSignalMode);
                end
            end
            
            % preparation des axes
            for kCheck = 1:length(Checkboxs2)
                cb = Checkboxs2(kCheck);
                if cb.Value
                    FiguresCheckboxs2(kCheck) = figure;
                    if multipleAxesDisplay
                        for kPlot = 1:nbPlots
                            axesFiguresCheckboxs2(kCheck, kPlot) =...
                                subplot0(nbPlots, 1, kPlot, axes(FiguresCheckboxs2(kCheck)));
                        end
                    else
                        axesFiguresCheckboxs2(kCheck, 1) = axes(FiguresCheckboxs2(kCheck));
                        for kPlot = 2:nbPlots
                            axesFiguresCheckboxs2(kCheck, kPlot) =...
                                axesFiguresCheckboxs2(kCheck, 1);
                        end
                        hold(axesFiguresCheckboxs2(kCheck, 1), 'on');
                    end
                end
            end
            
            
            % plot des ridges
            for kPlot = 1:length(ridges)
                if autocorrelationMode
                    XLimRidgePlot = tRy([1 end]);
                elseif rdtMode
                    XLimRidgePlot = Xrdt([1 end]);
                else
                    XLimRidgePlot = XLim;
                end
                
                ridge = ridges{kPlot};
                
                if  checkboxTimeAmplPlot.Value && ~isempty(plotAxes) && all(isvalid(plotAxes)) % plot de l'amplitude directement sur l'axe
                    [~, newTimeAmplPlots] = RidgeQtyPlot2(ridge, 'time', 'val', 'EvaluationFunctionY', 'abs',...
                        'Axes', plotAxes(kPlot), 'Grid', 'auto', 'RenameAxes', false);
                    timeAmplPlots = [timeAmplPlots, newTimeAmplPlots];
                end
                
                if checkboxTimeAmpl.Value % plot de l'amplitude
                    axRidges = RidgeQtyPlot2(ridge, 'time', 'val', 'EvaluationFunctionY', 'abs',...
                        'ScaleX', get(xscaleTimeAmpl, 'String'), 'ScaleY', get(yscaleTimeAmpl, 'String'),...
                        'Axes', axesFiguresCheckboxs2(1, kPlot), 'RenameAxes', true, 'NameX', 'Time [s]',...
                        'NameY', ['Amplitude [', signalUnit, ']'], 'XLim', XLimRidgePlot, 'Threshold', RidgeMinModu);
                    XLimDisplayObj.addAxes(axRidges, false(size(axRidges)));
                end
                if checkboxTimeFreq.Value % plot de la frequence
                    axRidges = RidgeQtyPlot2(ridge, 'time', freqRidgeName,...
                        'ScaleX', get(xscaleTimeFreq, 'String'), 'ScaleY', get(yscaleTimeFreq, 'String'),...
                        'Axes', axesFiguresCheckboxs2(2, kPlot), 'RenameAxes', true, 'NameX', 'Time [s]',...
                        'NameY', 'Frequency [Hz]', 'XLim', XLimRidgePlot, 'Threshold', RidgeMinModu);
                    XLimDisplayObj.addAxes(axRidges, false(size(axRidges)));
                end
                if checkboxTimeDamp.Value % plot de l'amortissement
                    axRidges = RidgeQtyPlot2(ridge, 'time', dampingRidgeName,...
                        'ScaleX', get(xscaleTimeDamp, 'String'), 'ScaleY', get(yscaleTimeDamp, 'String'),...
                        'Axes', axesFiguresCheckboxs2(3, kPlot), 'RenameAxes', true, 'NameX', 'Time [s]',...
                        'NameY', 'Damping', 'XLim', XLimRidgePlot, 'Threshold', RidgeMinModu);
                    XLimDisplayObj.addAxes(axRidges, false(size(axRidges)));
                end
                if checkboxAmplFreq.Value % plot de l'amplitude en fonction de la frequence
                    RidgeQtyPlot2(ridge, 'val', freqRidgeName, 'EvaluationFunctionX', 'abs',...
                        'ScaleX', get(xscaleAmplFreq, 'String'), 'ScaleY', get(yscaleAmplFreq, 'String'),...
                        'Axes', axesFiguresCheckboxs2(4, kPlot), 'RenameAxes', true,...
                        'NameX', ['Amplitude [', signalUnit, ']'], 'NameY', 'Frequency [Hz]', 'Threshold', RidgeMinModu);
                end
                if checkboxAmplDamp.Value % plot de l'amortissement
                    RidgeQtyPlot2(ridge, 'val', dampingRidgeName, 'EvaluationFunctionX', 'abs',...
                        'ScaleX', get(xscaleAmplDamp, 'String'), 'ScaleY', get(yscaleAmplDamp, 'String'),...
                        'Axes', axesFiguresCheckboxs2(5, kPlot), 'RenameAxes', true,...
                        'NameX', ['Amplitude [', signalUnit, ']'], 'NameY', 'Damping', 'Threshold', RidgeMinModu);
                end
                if checkboxTimePhase.Value % plot de la phase
                    axRidges = RidgeQtyPlot2(ridge, 'time', phaseRidgeName,...
                        'ScaleX', get(xscaleTimePhase, 'String'), 'ScaleY', get(yscaleTimePhase, 'String'),...
                        'Axes', axesFiguresCheckboxs2(6, kPlot), 'RenameAxes', true,...
                        'NameX', 'Time [s]', 'NameY', 'Phase (rad)', 'XLim', XLimRidgePlot, 'Threshold', RidgeMinModu);
                    XLimDisplayObj.addAxes(axRidges, false(size(axRidges)));
                end
            end
        end
        
        %% calcul des deformees
        if any([Checkboxs3.Value]) || any([Checkboxs4.Value])
            if multiSignalMode || (autocorrelationMode && autocorrelationSVDMode)
                if multiSignalMode
                    if rdtMode
                        xMode = Xrdt;
                        yMode = Yrdt;
                    else
                        xMode = x(1,:);
                        yMode = y;
                    end
                    [timeShapes, freqsShapes, shapesShapes, amplitudesShapes] = ...
                        getModesSingleRidge(xMode, yMode, Q, fmin, fmax, NbFreq,...
                        'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu,...
                        'StopWhenIncreasing', StopRidgeWhenIncreasing,...
                        'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects, 'MaxSlopeRidge', slopeRidge,...
                        'MotherWavelet', MotherWavelet, 'XLimRidge', XLimRidge, 'ctRidge', ctRidge,...
                        'FrequencyScale', FrequencyScale, 'MultiSignalModeAbsValue', multiSignalModeAbsValue);
                elseif autocorrelationMode
                    [timeShapesSVD, freqsShapesSVD, shapesShapesSVD, amplitudesShapesSVD] = ...
                        getModesCrossCorr(tRy, SVry, SVvectry, Q, fmin, fmax, NbFreq, autocorrelationNsvd,...
                        'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu,...
                        'StopWhenIncreasing', StopRidgeWhenIncreasing,...
                        'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects, 'MaxSlopeRidge', slopeRidge,...
                        'MotherWavelet', MotherWavelet, 'XLimRidge', XLimRidge, 'ctRidge', ctRidge,...
                        'FrequencyScale', FrequencyScale);
                end
                
                if autocorrelationMode
                    Nsvd = autocorrelationNsvd;
                else
                    Nsvd = 1;
                end
                
                for ksvd = 1:Nsvd
                    if autocorrelationMode
                        timeShapes = timeShapesSVD{ksvd};
                        freqsShapes = freqsShapesSVD{ksvd};
                        shapesShapes = shapesShapesSVD{ksvd};
                        amplitudesShapes = amplitudesShapesSVD{ksvd};
                    end
                    
                    
                    absAmplitudesShapes = {};
                    for kridge = 1:length(amplitudesShapes)
                        absAmplitudesShapes{end+1} = abs(amplitudesShapes{kridge});
                    end
                    
                    % axe x
                    XquantName = xAxisShapesnames{get(xAxisShapes, 'Value')};
                    if isequal(XquantName, 'time')
                        XquantLabel = 'Time [s]';
                        XquantScale = 'lin';
                        Xquantity = timeShapes;
                    elseif isequal(XquantName, 'ampl')
                        XquantLabel = ['Amplitude [', signalUnit, ']'];
                        XquantScale = 'lin';
                        Xquantity = absAmplitudesShapes;
                    elseif isequal(XquantName, 'log(ampl)')
                        XquantLabel = ['Amplitude [', signalUnit, ']'];
                        XquantScale = 'log';
                        Xquantity = absAmplitudesShapes;
                    elseif isequal(XquantName, 'freq')
                        XquantLabel = 'Frequence [Hz]';
                        XquantScale = 'lin';
                        Xquantity = freqsShapes;
                    end
                    
                    % classment ordre modes
                    meanFreqsModes = nan(1, length(timeShapes));
                    for kridge = 1:length(timeShapes)
                        % moyenne
                        if get(weightOptionShapesMean, 'Value')
                            weights = abs(amplitudesShapes{kridge});
                        else
                            weights = ones(size(amplitudesShapes{kridge}));
                        end
                        weights = weights / sum(weights);
                        meanFreqsModes(kridge) = sum(freqsShapes{kridge} .* weights);
                    end
                    
                    if false % tri des freqs
                        [meanFreqsModes, modesOrder] = sort(meanFreqsModes);
                    else
                        modesOrder = 1:length(timeShapes);
                    end
                    
                    % plots
                    for kmode = 1:length(modesOrder)
                        kridge = modesOrder(kmode);
                        figuresName = ['mode ', num2str(meanFreqsModes(kmode)), 'Hz (', plotAxesName, ')'];
                        if autocorrelationMode
                            figuresName = ['SVD', num2str(ksvd), ' ', figuresName];
                        elseif rdtMode
                            figuresName = ['RDT ', figuresName];
                        end
                        if checkboxRealShapes.Value
                            figShape = figure('Name', figuresName);
                            ax = axes(figShape);
                            plot(ax, Xquantity{kridge}, real(shapesShapes{kridge}));
                            xlabel(ax, XquantLabel);
                            ylabel(ax, 'Real part');
                            set(ax, 'XScale', XquantScale);
                        end
                        if checkboxImagShapes.Value
                            figShape = figure('Name', figuresName);
                            ax = axes(figShape);
                            plot(ax, Xquantity{kridge}, imag(shapesShapes{kridge}));
                            xlabel(ax, XquantLabel);
                            ylabel(ax, 'Imaginary part');
                            set(ax, 'XScale', XquantScale);
                        end
                        if checkboxModuleShapes.Value
                            figShape = figure('Name', figuresName);
                            ax = axes(figShape);
                            plot(ax, Xquantity{kridge}, abs(shapesShapes{kridge}));
                            xlabel(ax, XquantLabel);
                            ylabel(ax, 'Module');
                            set(ax, 'XScale', XquantScale);
                        end
                        if checkboxPhaseShapes.Value
                            figShape = figure('Name', figuresName);
                            ax = axes(figShape);
                            plot(ax, Xquantity{kridge}, angle(shapesShapes{kridge}));
                            xlabel(ax, XquantLabel);
                            ylabel(ax, 'Phase');
                            set(ax, 'XScale', XquantScale);
                        end
                        if checkboxAmplShapes.Value
                            figShape = figure('Name', figuresName);
                            ax = axes(figShape);
                            plot(ax, Xquantity{kridge}, abs(amplitudesShapes{kridge}));
                            xlabel(ax, XquantLabel);
                            ylabel(ax, 'Amplitude');
                            set(ax, 'XScale', XquantScale);
                        end
                        
                        % moyenne
                        if get(weightOptionShapesMean, 'Value')
                            weights = abs(amplitudesShapes{kridge});
                        else
                            weights = ones(size(amplitudesShapes{kridge}));
                        end
                        weights = weights / sum(weights);
                        meanShape = sum(shapesShapes{kridge} .* weights, 2);
                        meanFreq = sum(freqsShapes{kridge} .* weights);
                        
                        if checkboxRealShapesMean.Value
                            RealShapePlot(real(meanShape), figuresName)
                        end
                        if checkboxComplexShapesMean.Value
                            ComplexShapePlot(meanShape, figuresName)
                        end
                        if checkboxDispShapesMean.Value
                            disp(['mean shape, ', figuresName]);
                            disp(meanShape);
                        end
                        if checkboxDispFreqsMean.Value
                            disp(['mean frequency, ', figuresName]);
                            disp(meanFreq);
                        end
                        if checkboxAmplRegMean.Value
                            figShape = figure('Name', figuresName);
                            ax = axes(figShape);
                            plot(ax, timeShapes{kridge}, log(abs(amplitudesShapes{kridge})));
                            xlabel(ax, 'Time [s]');
                            ylabel(ax, 'Amplitude log');
                            hold(ax, 'on');
                            % linear regression
                            regResults = [ones(size(timeShapes{kridge})); timeShapes{kridge}].' \...
                                log(abs(amplitudesShapes{kridge})).';
                            plot(ax, timeShapes{kridge}, regResults(1) + regResults(2)*timeShapes{kridge},...
                                'r--');
                            % modal computations
                            decayRate = -regResults(2);
                            naturalFreq = sqrt(meanFreq^2 + (decayRate/(2*pi))^2);
                            dampingRatio = decayRate/(2*pi)/naturalFreq;
                            fprintf('damped frequency: %.10f Hz\n', meanFreq);
                            fprintf('decay rate: %.10f Hz\n', decayRate);
                            fprintf('-> natural frequency: %.10f Hz\n', naturalFreq);
                            fprintf('-> damping ratio: %.10f %%\n\n', 100*dampingRatio);
                        end
                    end
                end
            else
                warning('! multi signal or autocorrelation SVD modes only !');
            end
        end
        
        % changement du curseur
        set(fig, 'pointer', 'arrow');
        drawnow;
        
    end

%%

    function plotCrossCorr(Dx, y, NmaxLagCorr, forcePlot)
        % plot iif Ry change | forcePlot
        % forcePlot == 1 for cross-corr, == 2 for autocorr
        if nargin < 4
            forcePlot = isempty(tRy);
        end
        
        if isempty(tRy)
            tRy = Dx * (0:NmaxLagCorr);
            Ry = crossCorrelation(y, NmaxLagCorr, autocorrelationBias);
        end
        
        if forcePlot
            figRy = figure('Name', ['auto/cross-corr (', plotAxesName, ')']);
            axRy = axes(figRy);
            xlabel(axRy, 'Time [s]');
            if forcePlot == 1
                ylabel(axRy, ['Cross-correlation [', squaredSignalUnit, ']']);
            elseif forcePlot == 2
                ylabel(axRy, ['Autocorrelation [', squaredSignalUnit, ']']);
            end
            XLimDisplayObj.addAxesXLimRidges(axRy, true)
            hold (axRy, 'on');
            
            Ndof = size(y, 1);
            legendRij = {};
            for i = 1:Ndof
                plot(axRy, tRy, reshape(Ry(i, i, :), [1, size(Ry, 3)]));
                legendRij{end+1} = ['R', num2str(i), num2str(i)];
            end
            if forcePlot == 1
                for i = 1:Ndof
                    for j = 1:Ndof
                        if i ~= j
                            plot(axRy, tRy, reshape(Ry(i, j, :), [1, size(Ry, 3)]), ':');
                            legendRij{end+1} = ['R', num2str(i), num2str(j)];
                        end
                    end
                end
            end
            
%             xlim(axRy, [tRy(1), tRy(end)]);
            legend(axRy, legendRij{:});
            drawnow;
        end
    end


    function resetCrossCorr()
        tRy = [];
        Ry = [];
        SVry = [];
        SVvectry = [];
    end


%%

    function deletePlots()
        try
            delete(timeAmplPlots);
        catch
        end
        timeAmplPlots = [];
    end


%% hilbert

    function hilbertTransform()
        getXY();
        
        hilbertPlotAxes = plotAxes;
        if isempty(hilbertPlotAxes) || ~all(isvalid(hilbertPlotAxes))
            fhilb = figure;
            hilbertPlotAxes = [];
            for kPlot = 1:nbPlots
                hilbertPlotAxes(kPlot) = subplot0(nbPlots, 1, kPlot, axes(fhilb));
                if ~multipleAxesDisplay
                    hold(axes(hilbertPlotAxes(kPlot)), 'on');
                end
            end
        end
        
        for kPlot = 1:nbPlots
            hilb = hilbert(y(kPlot,:));
            hold(hilbertPlotAxes(kPlot), 'on');
            newPlot = plot(x(kPlot,:), abs(hilb), 'Parent', hilbertPlotAxes(kPlot));
            hold(hilbertPlotAxes(kPlot), 'off');
            timeAmplPlots = [timeAmplPlots, newPlot];
        end
            
    end


%% fourier

    function fourierTransformDisplay()
        getXY();
        
        fminNew = eval(get(editfmin, 'String'));
        fmaxNew = eval(get(editfmax, 'String'));
        if fminNew ~= fmin || fmaxNew ~= fmax 
            resetCrossCorr();
        end
        fmin = fminNew;
        fmax = fmaxNew;
        
        if autocorrelationMode
            Dx = (x(1, end) - x(1, 1))/(size(x, 2)-1);
            NmaxLagCorr = floor(maxLagCorr/Dx);
            plotCrossCorr(Dx, y, NmaxLagCorr);
        end
        
        if ~autocorrelationMode || (autocorrelationMode && ~autocorrelationFourierSVDMode)
            ffourier = figure;
            fourierPlotAxes = [];
            
            nbAxes = nbPlots;
            
            if multipleAxesDisplay %création des axes où sont plot les courbes
                for kPlot = 1:nbAxes
                    fourierPlotAxes(kPlot) = subplot0(nbAxes, 1, kPlot, axes(ffourier));
                    set(fourierPlotAxes(kPlot), 'XScale', FrequencyScale, 'YScale', 'lin');
                    set(fourierPlotAxes(kPlot), 'Xlim', [fmin fmax]);
                end
            else
                axesffourier = axes(ffourier);
                hold(axesffourier, 'on');
                set(axesffourier, 'XScale', FrequencyScale, 'YScale', 'lin');
                for kPlot = 1:nbAxes
                    fourierPlotAxes(kPlot) = axesffourier;
                    set(fourierPlotAxes(kPlot), 'Xlim', [fmin fmax]);
                end
            end
            
            if ~multiSignalMode % plot des courbes
                for kPlot = 1:nbAxes
                    if autocorrelationMode && ~autocorrelationFourierSVDMode
                        Xfour = tRy;
                        Yfour = Ry(kPlot,kPlot,:); % fft cross-corr TODO
                        Yfour = transpose(Yfour(:));
                    elseif rdtMode
                        Xfour = Xrdt;
                        Yfour = Yrdt(kPlot,:);
                    else
                        Xfour = x(kPlot,:);
                        Yfour = y(kPlot,:);
                    end
                    Yfour = [Yfour, zeros(1, ZeroPaddingFourier*length(Yfour))];
                    Xfour = mean(diff(Xfour)) * (0:length(Yfour)-1);
                    Tfour = mean(diff(Xfour)) * length(Xfour);
                    
                    [freqs, four] = fourierTransform(Xfour, Yfour,...
                        'Averaging', FourierAveraging, 'AveragingNb', FourierAveragingNb,...
                        'Window', FourierWindow, 'WindowParams', FourierWindowParams);
                    
                    hold(fourierPlotAxes(kPlot), 'on');
                    if isequal(FourierScale, 'lin')
                        pltFourier = plot(fourierPlotAxes(kPlot), freqs, abs(four));
                    elseif isequal(FourierScale, 'squared')
                        pltFourier = plot(fourierPlotAxes(kPlot), freqs, abs(four).^2);
                    elseif isequal(FourierScale, 'log')
                        pltFourier = plot(fourierPlotAxes(kPlot), freqs, abs(four));
                        set(fourierPlotAxes(kPlot), 'YScale', 'log');
                    elseif isequal(FourierScale, 'spectral density (lin)')
                        pltFourier = plot(fourierPlotAxes(kPlot), freqs, abs(four).^2 * Tfour);
                    elseif isequal(FourierScale, 'spectral density (log)')
                        pltFourier = plot(fourierPlotAxes(kPlot), freqs, abs(four).^2 * Tfour);
                        set(fourierPlotAxes(kPlot), 'YScale', 'log');
                    elseif isequal(FourierScale, 'phase')
                        pltFourier = plot(fourierPlotAxes(kPlot), freqs, angle(four));
                    end
                    hold(fourierPlotAxes(kPlot), 'off');
                end
            else % muli signal mode
                FourierTot = 0;
                for kPlot = 1:nbAxes
                    if autocorrelationMode && ~autocorrelationFourierSVDMode
                        Xfour = tRy;
                        Yfour = Ry(kPlot,kPlot,:);
                        Yfour = transpose(Yfour(:));
                    elseif rdtMode
                        Xfour = Xrdt;
                        Yfour = Yrdt(kPlot,:);
                    else
                        Xfour = x(kPlot,:);
                        Yfour = y(kPlot,:);
                    end
                    Yfour = [Yfour, zeros(1, ZeroPaddingFourier*length(Yfour))];
                    Xfour = mean(diff(Xfour)) * (0:length(Yfour)-1);
                    
                    [freqs, four] = fourierTransform(Xfour, Yfour,...
                        'Averaging', FourierAveraging, 'AveragingNb', FourierAveragingNb,...
                        'Window', FourierWindow, 'WindowParams', FourierWindowParams);
                    
                    if ~multiSignalModeAbsValue
                        FourierTot = FourierTot + four.^2;
                    else
                        FourierTot = FourierTot + abs(four).^2;
                    end
                end
                four = sqrt(FourierTot);
                hold(fourierPlotAxes(kPlot), 'on');
                if isequal(FourierScale, 'lin')
                    plot(fourierPlotAxes(kPlot), freqs, abs(four));
                elseif isequal(FourierScale, 'squared')
                    plot(fourierPlotAxes(kPlot), freqs, abs(four).^2);
                elseif isequal(FourierScale, 'log')
                    plot(fourierPlotAxes(kPlot), freqs, abs(four));
                    set(fourierPlotAxes(kPlot), 'YScale', 'log');
                elseif isequal(FourierScale, 'spectral density (lin)')
                    plot(fourierPlotAxes(kPlot), freqs, abs(four).^2 * Tfour);
                elseif isequal(FourierScale, 'spectral density (log)')
                    plot(fourierPlotAxes(kPlot), freqs, abs(four).^2 * Tfour);
                    set(fourierPlotAxes(kPlot), 'YScale', 'log');
                elseif isequal(FourierScale, 'phase')
                    plot(fourierPlotAxes(kPlot), freqs, angle(four));
                end
                hold(fourierPlotAxes(kPlot), 'off');
            end
            
        else % autocorr
            if autocorrelationFourierSVDMode
                [freqs, SVfftrx] = svdFFT(tRy, Ry, autocorrelationNsvd,...
                    'Window', FourierWindow, 'WindowParams', FourierWindowParams);
                
                fourierFig = figure;
                fourierPlotAxes = axes(fourierFig);
                hold(fourierPlotAxes, 'on');
                for ksv = 1:autocorrelationNsvd
                    fftRx = SVfftrx{ksv};
                    Tfour = length(tRy) * mean(diff(tRy));
                    
                    if isequal(FourierScale, 'lin')
                        pltSVDfourier = plot(fourierPlotAxes, freqs, abs(fftRx));
                    elseif isequal(FourierScale, 'squared')
                        pltSVDfourier = plot(fourierPlotAxes, freqs, abs(fftRx).^2);
                    elseif isequal(FourierScale, 'log')
                        pltSVDfourier = plot(fourierPlotAxes, freqs, abs(fftRx));
                        set(fourierPlotAxes, 'YScale', 'log');
                    elseif isequal(FourierScale, 'spectral density (lin)')
                        pltSVDfourier = plot(fourierPlotAxes, freqs, abs(fftRx).^2 * Tfour);
                    elseif isequal(FourierScale, 'spectral density (log)')
                        pltSVDfourier = plot(fourierPlotAxes, freqs, abs(fftRx).^2 * Tfour);
                        set(fourierPlotAxes, 'YScale', 'log');
                    elseif isequal(FourierScale, 'phase')
                        pltSVDfourier = plot(fourierPlotAxes, freqs, angle(fftRx));
                    end
                    
                    set(pltSVDfourier, 'DisplayName', sprintf('SV%u', ksv));
                    
                    set(fourierPlotAxes, 'Xlim', [fmin fmax]);
                    set(fourierPlotAxes, 'XScale', FrequencyScale);
                end
                legend(fourierPlotAxes);
            else % deprecated
                fourierFig = figure;
                fourierPlotAxes = axes(fourierFig);
                hold(fourierPlotAxes, 'on');
                
                Tfour = mean(diff(tRy)) * length(tRy);
                
                legendRij = cell(1, nbPlots^2);
                for i = 1:nbPlots
                    for j = 1:nbPlots
                        Yfour = reshape(Ry(i, j, :), [1, size(Ry, 3)]);
                        Yfour = [Yfour, zeros(1, ZeroPaddingFourier*length(Yfour))];
                        Xfour = tRy;
                        Xfour = mean(diff(Xfour)) * (0:length(Xfour)-1);
                        
                        [freqs, fftRx] = fourierTransform(Xfour, Yfour,...
                            'Averaging', FourierAveraging, 'AveragingNb', FourierAveragingNb,...
                        'Window', FourierWindow, 'WindowParams', FourierWindowParams);
                        
                        if isequal(FourierScale, 'lin')
                            fourierPlot = plot(fourierPlotAxes, freqs, abs(fftRx));
                        elseif isequal(FourierScale, 'squared')
                            fourierPlot = plot(fourierPlotAxes, freqs, abs(fftRx).^2);
                        elseif isequal(FourierScale, 'log')
                            fourierPlot = plot(fourierPlotAxes, freqs, abs(fftRx));
                            set(fourierPlotAxes, 'YScale', 'log');
                        elseif isequal(FourierScale, 'spectral density (lin)')
                            fourierPlot = plot(fourierPlotAxes, freqs, abs(fftRx).^2 * Tfour);
                        elseif isequal(FourierScale, 'spectral density (log)')
                            fourierPlot = plot(fourierPlotAxes, freqs, abs(fftRx).^2 * Tfour);
                            set(fourierPlotAxes, 'YScale', 'log');
                        elseif isequal(FourierScale, 'phase')
                            fourierPlot = plot(fourierPlotAxes, freqs, angle(fftRx));
                        end
                        
                        if i ~= j
                            set(fourierPlot, 'LineStyle', ':');
                        end
                        
                        legendRij{(i-1)*nbPlots + j} = ['R', num2str(i), num2str(j)];
                    end
                end
                
                set(fourierPlotAxes, 'Xlim', [fmin fmax]);
                set(fourierPlotAxes, 'XScale', FrequencyScale);
                legend(fourierPlotAxes, legendRij);
                
                hold(fourierPlotAxes, 'off');
                
            end
        end
        
        for kPlot = 1:length(fourierPlotAxes)
            xlabel(fourierPlotAxes(kPlot), 'Frequency [Hz]');
            
            if autocorrelationMode
                signalUnitFourier = squaredSignalUnit;
                squaredSignalUnitFourier = ['(', squaredSignalUnit, ')²'];
            else
                signalUnitFourier = signalUnit;
                squaredSignalUnitFourier = squaredSignalUnit;
            end
            
            switch FourierScale
                case {'lin', 'log'}
                    ylabel(fourierPlotAxes(kPlot), ['Amplitude [', signalUnitFourier, ']']);
                case 'squared'
                    ylabel(fourierPlotAxes(kPlot), ['Squared amplitude [', squaredSignalUnitFourier, ']']);
                case {'spectral density (lin)', 'spectral density (log)'}
                    ylabel(fourierPlotAxes(kPlot), ['Spectral density [', squaredSignalUnitFourier, '/Hz]']);
                case 'phase'
                    ylabel(fourierPlotAxes(kPlot), 'Phase [rad]');
            end
        end
    end


buttonWavelet.Callback = @(~,~) show();

deleteButton.Callback = @(~,~) deletePlots();

buttonHilbert.Callback = @(~,~) hilbertTransform();

buttonFourier.Callback = @(~,~) fourierTransformDisplay();


end