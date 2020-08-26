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
defaultCtEdgeEffects = 3;
defaultZeroPaddingFourier = 0;
defaultMultiSignalMode = false;
defaultAutocorrelationMode = false;
defaultMaxLagCorr = nan;
defaultAutocorrelationNsvd = 1;
defaultWvltScale = 'log';
defaultFourierScale = 'lin';
defaultXLim = nan;
defaultWvltAxesTitle = '';
defaultComplexShapePlot = @complexShapePlot1;
defaultRealShapePlot = @realShapePlot1;
defaultMotherWavelet = 'cauchy';



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
addParameter(p,'CtEdgeEffects', defaultCtEdgeEffects);
addParameter(p,'ZeroPaddingFourier', defaultZeroPaddingFourier);
addParameter(p,'MultiSignalMode', defaultMultiSignalMode);
addParameter(p,'AutocorrelationMode', defaultAutocorrelationMode);
addParameter(p,'AutocorrelationMaxLag', defaultMaxLagCorr);
addParameter(p,'AutocorrelationNsvd', defaultAutocorrelationNsvd);
addParameter(p,'WvltScale', defaultWvltScale);
addParameter(p,'FourierScale', defaultFourierScale);
addParameter(p,'XLim', defaultXLim);
addParameter(p, 'WvltAxesTitle', defaultWvltAxesTitle);
addParameter(p, 'ComplexShapePlot', defaultComplexShapePlot);
addParameter(p, 'RealShapePlot', defaultRealShapePlot);
addParameter(p, 'MotherWavelet', defaultMotherWavelet);

parse(p, varargin{:})

fig = p.Results.Parent;
if fig == 0
    fig = figure('Name', 'Wavelet Menu', 'numbertitle', 'off');
    fig.Units = 'characters';
    fig.Position(3) = 110;
    fig.Position(4) = 26;
    fig.MenuBar = 'none';
%     fig.ToolBar = 'none';
end

    function M = cellmat2mat(M)
        if isa(M, 'cell')
            for iM = 1:length(M)
                M{iM} = M{iM}';
            end
            M = [M{:}]';
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
    plotAxes = 0;
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
nbPlots = size(getX(), 1);


fmin = p.Results.fmin;
fmax = p.Results.fmax;
NbFreq = p.Results.NbFreq;
Q = p.Results.Q;
MaxRidges = p.Results.MaxRidges;
MaxParallelRidges = p.Results.MaxParallelRidges;
MaxSlopeRidge = p.Results.MaxSlopeRidge;
multipleAxesDisplay = p.Results.MultipleAxesDisplay;
RidgeMinModu = p.Results.RidgeMinModu;
ctEdgeEffects = p.Results.CtEdgeEffects;
ZeroPaddingFourier = p.Results.ZeroPaddingFourier;

multiSignalMode = p.Results.MultiSignalMode;
autocorrelationMode = p.Results.AutocorrelationMode;
tRy = nan; Ry = nan; SVry = nan; SVvectry = nan;
resetCrossCorr();

autocorrelationNsvd = p.Results.AutocorrelationNsvd;
maxLagCorr = p.Results.AutocorrelationMaxLag;
WvltScale = p.Results.WvltScale;
FourierScale = p.Results.FourierScale;
XLim = p.Results.XLim;
wvltAxesTitle = p.Results.WvltAxesTitle;
ComplexShapePlot = p.Results.ComplexShapePlot;
RealShapePlot = p.Results.RealShapePlot;
MotherWavelet = p.Results.MotherWavelet;

x0 = getX();
if isnan(XLim)
    Xmin = x0(1);
    Xmax = x0(end);
else
    Xmin = XLim(1);
    Xmax = XLim(2);
end

if isnan(maxLagCorr)
    maxLagCorr = Xmax-Xmin;
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


%% bouton ondelettes et panneaux param et sorties

%ondelette
waveletPan = uipanel('Parent',fig, 'Units', 'normalized');
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

strPR = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
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

checkboxPhase = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'phase plot', 'Value', false);



Checkboxs1 = [checkboxGeneral, checkboxModule, checkboxPhase];

if ~isequal(plotAxes, 0)
    checkboxTimeAmplPlot = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
        'String', 'time, ampl on plot', 'Value', false);
    Checkboxs1 = [Checkboxs1, checkboxTimeAmplPlot];
    timeAmplPlots = [];
    deleteButton = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','pushbutton',...
        'String', 'delete');
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

checkboxAmplFreq = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'ampl, freq', 'Value', false);
xscaleAmplFreq = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleAmplFreq = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

checkboxTimeDamp = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'time, damping', 'Value', false);
xscaleTimeDamp = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleTimeDamp = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
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

Checkboxs2 = [checkboxTimeAmpl, checkboxTimeFreq, checkboxAmplFreq, checkboxTimeDamp, checkboxAmplDamp,...
    checkboxTimePhase];
XScales = [xscaleTimeAmpl, xscaleTimeFreq, xscaleAmplFreq, xscaleTimeDamp, xscaleAmplDamp, xscaleTimePhase];
YScales = [yscaleTimeAmpl, yscaleTimeFreq, yscaleAmplFreq, yscaleTimeDamp, yscaleAmplDamp, yscaleTimePhase];
n1 = length(Checkboxs1);
n2 = length(Checkboxs2);
n = n1+n2;
for k=1:n1
    Checkboxs1(k).Position = [0.01, 0.01+(n-k)/n, 0.48, 1/n-0.02];
end
for k=(n1+1:n1+n2)
    Checkboxs2(k-n1).Position = [0.01, 0.01+(n-k)/n, 0.48, 1/n-0.02];
end
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

LogScales = [yscaleTimeAmpl, xscaleAmplFreq];
for scale = LogScales
    set(scale, 'Value', 1);
    set(scale, 'String', scalesNames{scale.Value+1});
end





%% sorties mode shapes

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
Checkboxs3 = [checkboxRealShapes, checkboxImagShapes, checkboxModuleShapes, checkboxPhaseShapes];



strShapesMean = uicontrol('Parent', shapesPan, 'Units', 'normalized', 'Style', 'text',...
    'String', 'average:', 'HorizontalAlignment', 'left');
weightOptionShapesMean = uicontrol('Parent', shapesPan, 'Units', 'normalized', 'Style','togglebutton',...
    'String', 'weighted (ampl)', 'Value', false);

checkboxRealShapesMean = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'plot real', 'Value', false);
checkboxComplexShapesMean = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'plot complex', 'Value', false);
checkboxDispShapesMean = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'disp complex', 'Value', false);
checkboxDispFreqsMean = uicontrol('Parent',shapesPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'frequency', 'Value', false);
Checkboxs4 = [checkboxRealShapesMean, checkboxComplexShapesMean, checkboxDispShapesMean, checkboxDispFreqsMean];



n3 = ceil(length(Checkboxs3)/2);
n4 = ceil(length(Checkboxs4)/2);
n = n3 + n4 + 5.5;

strxAxisShapes.Position = [0.05, (n-1.5)/n, 0.94, 1/n-0.02];
xAxisShapes.Position = [0.01, (n-2.5)/n+0.01, 0.48, 1/n-0.02];

stryAxisShapes.Position = [0.05, (n-4)/n, 0.94, 1/n-0.02];
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
MotherWaveletNames = {'cauchy', 'morlet', 'harmonic', 'littlewood-paley'};

motherWaveletMenu = uimenu(paramMenu, 'Text','Mother Wavelet');
motherWaveletMenuChoices(1) = uimenu(motherWaveletMenu, 'Text', 'Cauchy', 'Checked' ,'on');
motherWaveletMenuChoices(2) = uimenu(motherWaveletMenu, 'Text', 'Morlet');
motherWaveletMenuChoices(3) = uimenu(motherWaveletMenu, 'Text', 'Harmonic');
motherWaveletMenuChoices(4) = uimenu(motherWaveletMenu, 'Text', 'Littlewood-Paley');
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

%autres valeurs par défault
freqRidgeName = 'freq';
phaseRidgeName = 'pha2';
dampingRidgeName = 'damping3';

freqRidgeNames = {'freq', 'freq2'};
phaseRidgeNames = {'pha', 'pha2'};
dampingRidgeNames = {'damping', 'damping2', 'damping3'};
WvltScaleNames = {'lin', 'log10'};
FourierScaleNames = {'lin', 'squared', 'log', 'phase'};

%freq
freqMenu = uimenu(paramMenu, 'Text','Frequency');
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
phaseMenu = uimenu(paramMenu, 'Text','Phase');
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
dampingMenu = uimenu(paramMenu, 'Text','Damping');
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

% wavelet scale
WvltScaleMenu = uimenu(paramMenu, 'Text','Wavelet Scale');
WvltScaleMenuChoices(1) = uimenu(WvltScaleMenu, 'Text', 'lin');
WvltScaleMenuChoices(2) = uimenu(WvltScaleMenu, 'Text', 'log', 'Checked' ,'on');
    function selectWvltScaleMenu(kchoice)
        for kchoices = 1:length(WvltScaleMenuChoices)
            set(WvltScaleMenuChoices(kchoices), 'Checked', 'off');
        end
        set(WvltScaleMenuChoices(kchoice), 'Checked', 'on');
        
        WvltScale = WvltScaleNames{kchoice};
    end
set(WvltScaleMenuChoices(1), 'CallBack', @(~,~) selectWvltScaleMenu(1));
set(WvltScaleMenuChoices(2), 'CallBack', @(~,~) selectWvltScaleMenu(2));

% fourier scale
FourierScaleMenu = uimenu(paramMenu, 'Text','Fourier Scale');
FourierScaleMenuChoices(1) = uimenu(FourierScaleMenu, 'Text', 'lin', 'Checked', 'on');
FourierScaleMenuChoices(2) = uimenu(FourierScaleMenu, 'Text', 'squared');
FourierScaleMenuChoices(3) = uimenu(FourierScaleMenu, 'Text', 'log');
FourierScaleMenuChoices(4) = uimenu(FourierScaleMenu, 'Text', 'phase');
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

% Xlim
XlimMenu = uimenu(paramMenu, 'Text','Set Xlim');
    function setXlim()
        prompt = {'Enter Xmin :', 'Enter Xmax :'};
        dlgtitle = 'Input Xlim';
        dims = [1 35];
        definput = {num2str(Xmin), num2str(Xmax)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        try
            Xmin = str2double(answer{1});
            Xmax = str2double(answer{2});
        catch
        end
    end
set(XlimMenu, 'CallBack', @(~,~) setXlim);

% %zero padding fourier
% zeroPaddingFourierMenu = uimenu(paramMenu, 'Text','Zero padding Fourier');
%     function setZeroPaddingFourier()
%         prompt = {'Enter zero padding Fourier :'};
%         dlgtitle = 'Input zero padding Fourier';
%         dims = [1 35];
%         definput = {num2str(ZeroPaddingFourier)};
%         answer = inputdlg(prompt,dlgtitle,dims,definput);
%         try
%             ZeroPaddingFourier = str2double(answer{1});
%             ZeroPaddingFourier = max(round(ZeroPaddingFourier), 0);
%         catch
%         end
%     end
% set(zeroPaddingFourierMenu, 'CallBack', @(~,~) setZeroPaddingFourier);

%multipleAxesDisplay

multipleAxesDisplayMenu = uimenu(paramMenu, 'Text','Multiple axes', 'Checked', multipleAxesDisplay);
    function switchMultipleAxesDisplay(status)
        multipleAxesDisplayMenu.Checked = status;
        setMultipleAxesDisplay(status);
        if status
            switchMultiSignalModeDisplay(false);
            switchAutocorrelationModeDisplay(false)
        end
    end

multipleAxesDisplayMenu.MenuSelectedFcn = @(~, ~) switchMultipleAxesDisplay(~strcmp(multipleAxesDisplayMenu.Checked, 'on'));

%multiSignalMode

multiSignalModeMenu = uimenu(paramMenu, 'Text','Multi signal mode', 'Checked', multiSignalMode);
    function switchMultiSignalModeDisplay(status)
        multiSignalModeMenu.Checked = status;
        multiSignalMode = status;
        if status
            switchMultipleAxesDisplay(false);
            switchAutocorrelationModeDisplay(false);
        end
    end

multiSignalModeMenu.MenuSelectedFcn = @(~, ~) switchMultiSignalModeDisplay(~strcmp(multiSignalModeMenu.Checked, 'on'));

%autocorrelationMode

autocorrelationModeMenu = uimenu(paramMenu, 'Text','Cross-corr mode', 'Checked', autocorrelationMode);

autocorrelationParamsMenu = uimenu(paramMenu, 'Text', 'Cross-corr params');
autocorrelationDisplayMenu = uimenu(autocorrelationParamsMenu, 'Text', 'plot cross-corr');
autocorrelationMaxLagMenu = uimenu(autocorrelationParamsMenu, 'Text', 'set max lag');
autocorrelationFourierMenu = uimenu(autocorrelationParamsMenu, 'Text', 'fourier svd', 'Checked', 'on');

    function switchAutocorrelationModeDisplay(status)
        autocorrelationModeMenu.Checked = status;
        autocorrelationMode = status;
        set(autocorrelationParamsMenu, 'Enable', status);
        
        if status
            switchMultipleAxesDisplay(false);
            switchMultiSignalModeDisplay(false);
        end
    end
autocorrelationModeMenu.MenuSelectedFcn = @(~, ~) switchAutocorrelationModeDisplay(~strcmp(autocorrelationModeMenu.Checked, 'on'));

    function displayCrossCorr()
        getXY();
        Dx = (x(1, end) - x(1, 1))/(size(x, 2)-1);
        NmaxLagCorr = floor(maxLagCorr/Dx);
        plotCrossCorr(Dx, y, NmaxLagCorr, true);
    end
set(autocorrelationDisplayMenu, 'CallBack', @(~,~) displayCrossCorr);

    function setMaxLagCorr()
        answer = inputdlg({'Enter max lag :'}, 'Input MaxLag', [1 35], {num2str(maxLagCorr)});
        try
            maxLagCorr = str2double(answer{1});
        catch
            return
        end
        resetCrossCorr();
    end
set(autocorrelationMaxLagMenu, 'CallBack', @(~,~) setMaxLagCorr);

    function switchAutocorrelationFourierMenu()
        autocorrelationFourierMenu.Checked = ~strcmp(autocorrelationFourierMenu.Checked, 'on');
    end
set(autocorrelationFourierMenu, 'CallBack', @(~,~) switchAutocorrelationFourierMenu());

%set

switchMultiSignalModeDisplay(multiSignalMode);
switchAutocorrelationModeDisplay(autocorrelationMode);

%% menus regression & plot extract

toolsMenu = uimenu(fig, 'Text', 'Tools');

% filtering

filteringMenu = uimenu(toolsMenu,'Text','Filtering');

filteringMenu.MenuSelectedFcn = @(~, ~) LinearFilterMenu(WaveletPlot);

% regressions

regMenu = uimenu(toolsMenu,'Text','Regression');

regConstant = uimenu(regMenu, 'Text','Constant');
regExp = uimenu(regMenu, 'Text','Exp');

regConstant.MenuSelectedFcn = @(~, ~) RegressionMenu('Equation', 'c', 'Param', 'c', 'Param0', 0);
regExp.MenuSelectedFcn = @(~, ~) RegressionMenu('Equation', 'a*exp(-lambda*x)', 'Param', 'a  lambda',...
    'Param0', [1 1], 'Fit', 'log(y)');

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


%%

    function getXY()
        X = getX();
        Y = getY();
        
        for line = 1:size(X, 1)
            Y2(line, :) = Y(line, X(line, :)>=Xmin & X(line, :)<=Xmax);
            X2(line, :) = X(line, X(line, :)>=Xmin & X(line, :)<=Xmax);
        end
        
        %%%
        %test X identique
        %%%
        
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
        
        maxR = eval(get(editmaxR, 'String')); %nombre max de ridges
        PR = eval(get(editPR, 'String')); %nombre max de ridges parallèles
        slopeRidge = eval(get(editSL, 'String')); %max slope ridge
        getXY();
        Dx = (x(1, end) - x(1, 1))/(size(x, 2)-1);
        
        % cross corr
        if autocorrelationMode
            NmaxLagCorr = floor(maxLagCorr/Dx);
            plotCrossCorr(Dx, y, NmaxLagCorr);
            if isempty(SVry)
                [SVry, SVvectry] = svdCWT(tRy, Ry, fmin, fmax, NbFreq, Q, autocorrelationNsvd);
            end
        end
        
        %% plot de la transformee
        if checkboxGeneral.Value
            if ~multiSignalMode
                for kPlot = 1:nbPlots
                    WvltPlot(x(kPlot,:), y(kPlot,:), linspace(fmin,fmax,NbFreq), Q, 'ctEdgeEffects', ctEdgeEffects,...
                        'ctZeroPadding', ctEdgeEffects);
                end
            end
        end
        
        if checkboxModule.Value || checkboxPhase.Value
            if ~multiSignalMode && ~autocorrelationMode
                for kPlot = 1:nbPlots
                    wavelet = WvltComp(x(kPlot,:), y(kPlot,:), linspace(fmin,fmax,NbFreq), Q, 'ct', ctEdgeEffects, 'MotherWavelet', MotherWavelet);
                    if checkboxModule.Value
                        WvltPlot2(x(kPlot,:), linspace(fmin,fmax,NbFreq), wavelet, 'module', Q, ctEdgeEffects,...
                            WvltScale, ['Q=', num2str(Q),';scale:', WvltScale], wvltAxesTitle);
                    end
                    if checkboxPhase.Value
                        WvltPlot2(x(kPlot,:), linspace(fmin,fmax,NbFreq), wavelet, 'phase', Q, ctEdgeEffects,...
                            WvltScale, ['Q=', num2str(Q),';scale:', WvltScale], wvltAxesTitle);
                    end
                end
                
            elseif multiSignalMode
                wavelet = 0; % calcul de la somme des carrés de transformées
                for kPlot = 1:nbPlots
                    wavelet = wavelet +...
                        WvltComp(x(kPlot,:), y(kPlot,:), linspace(fmin,fmax,NbFreq), Q, 'ct', ctEdgeEffects, 'MotherWavelet', MotherWavelet).^2;
                end
                %wavelet = sqrt(wavelet);
                
                if checkboxModule.Value
                    WvltPlot2(x(kPlot,:), linspace(fmin,fmax,NbFreq), wavelet,...
                        'module', Q, ctEdgeEffects, WvltScale,...
                        ['sum_wvlt^2;Q=', num2str(Q),';scale:', WvltScale], wvltAxesTitle);
                end
                if checkboxPhase.Value
                    WvltPlot2(x(kPlot,:), linspace(fmin,fmax,NbFreq), wavelet,...
                        'phase', Q, ctEdgeEffects, WvltScale,...
                        ['sum_wvlt^2;Q=', num2str(Q),';scale:', WvltScale], wvltAxesTitle);
                end
                
            elseif autocorrelationMode
                for ksvd = 1:autocorrelationNsvd
                    if checkboxModule.Value
                        WvltPlot2(tRy, linspace(fmin,fmax,NbFreq), SVry{ksvd},...
                            'module', Q, ctEdgeEffects, WvltScale,...
                            ['xcorr->CWT->SVD;Q=', num2str(Q),';scale:', WvltScale], wvltAxesTitle);
                    end
                    if checkboxPhase.Value
                        WvltPlot2(tRy, linspace(fmin,fmax,NbFreq), SVry{ksvd},...
                            'phase', Q, ctEdgeEffects, WvltScale,...
                            ['sum_wvlt^2;Q=', num2str(Q),';scale:', WvltScale], wvltAxesTitle);
                    end
                end
            end
        end
        
        %% calcul des ridges
        if any([Checkboxs2.Value]) || (exist('checkboxTimeAmplPlot', 'var') && checkboxTimeAmplPlot.Value)
            
            if ~multiSignalMode && ~autocorrelationMode
                ridges = cell(1, nbPlots);
                for kPlot = 1:nbPlots
                    ridges{kPlot} = RidgeExtract(x(kPlot,:), y(kPlot,:), Q, fmin, fmax, NbFreq,...
                        'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu,...
                        'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects, 'MaxSlopeRidge', slopeRidge,...
                        'MotherWavelet', MotherWavelet);
                end
                
            elseif multiSignalMode
                wavelet = 0;
                for kPlot = 1:nbPlots
                    wavelet = wavelet +...
                        WvltComp(x(kPlot,:), y(kPlot,:), linspace(fmin,fmax,NbFreq), Q, 'ct', ctEdgeEffects, 'MotherWavelet', MotherWavelet).^2;
                end
                
                ridges = cell(1, 1);
                ridges{1} = RidgeExtract(x(kPlot,:), nan, Q, fmin, fmax, NbFreq,...
                    'Wavelet', wavelet,...
                    'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu^2,...
                    'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects, 'MaxSlopeRidge', slopeRidge,...
                    'MotherWavelet', MotherWavelet);
            
            elseif autocorrelationMode
                ridges = cell(1, autocorrelationNsvd);
                for ksvd = 1:autocorrelationNsvd
                    ridges{ksvd} = RidgeExtract(tRy, nan, Q, fmin, fmax, NbFreq, 'Wavelet', SVry{ksvd},...
                        'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu,...
                        'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects, 'MaxSlopeRidge', slopeRidge,...
                        'MotherWavelet', MotherWavelet);
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
                    XLimRidge = [tRy(1), tRy(end)];
                else
                    XLimRidge = [x(kPlot,1), x(kPlot,end)];
                end
                
                ridge = ridges{kPlot};
                
                if ~isequal(plotAxes, 0) && checkboxTimeAmplPlot.Value % plot de l'amplitude directement sur l'axe
                    newTimeAmplPlots = RidgeQtyPlot2(ridge, 'time', 'val', 'EvaluationFunctionY', 'abs',...
                        'Axes', plotAxes(kPlot), 'Grid', 'auto', 'RenameAxes', false, 'SquaredCWT', multiSignalMode);
                    timeAmplPlots = [timeAmplPlots, newTimeAmplPlots];
                end
                
                if checkboxTimeAmpl.Value % plot de l'amplitude
                    RidgeQtyPlot2(ridge, 'time', 'val', 'EvaluationFunctionY', 'abs',...
                        'ScaleX', get(xscaleTimeAmpl, 'String'), 'ScaleY', get(yscaleTimeAmpl, 'String'),...
                        'Axes', axesFiguresCheckboxs2(1, kPlot), 'RenameAxes', true, 'NameX', 'Time (s)', 'NameY', 'Amplitude (m/s^2)',...
                        'XLim', XLimRidge, 'SquaredCWT', multiSignalMode);
                end
                if checkboxTimeFreq.Value % plot de la frequence
                    RidgeQtyPlot2(ridge, 'time', freqRidgeName,...
                        'ScaleX', get(xscaleTimeFreq, 'String'), 'ScaleY', get(yscaleTimeFreq, 'String'),...
                        'Axes', axesFiguresCheckboxs2(2, kPlot), 'RenameAxes', true, 'NameX', 'Time (s)', 'NameY', 'Frequency (Hz)',...
                        'XLim', XLimRidge, 'SquaredCWT', multiSignalMode);
                end
                if checkboxAmplFreq.Value % plot de l'amplitude en fonction de la frequence
                    RidgeQtyPlot2(ridge, 'val', freqRidgeName, 'EvaluationFunctionX', 'abs',...
                        'ScaleX', get(xscaleAmplFreq, 'String'), 'ScaleY', get(yscaleAmplFreq, 'String'),...
                        'Axes', axesFiguresCheckboxs2(3, kPlot), 'RenameAxes', true, 'NameX', 'Amplitude (m/s^2)', 'NameY', 'Frequency (Hz)',...
						'SquaredCWT', multiSignalMode);
                end
                if checkboxTimeDamp.Value % plot de l'amortissement
                    RidgeQtyPlot2(ridge, 'time', dampingRidgeName,...
                        'ScaleX', get(xscaleTimeDamp, 'String'), 'ScaleY', get(yscaleTimeDamp, 'String'),...
                        'Axes', axesFiguresCheckboxs2(4, kPlot), 'RenameAxes', true, 'NameX', 'Time (s)', 'NameY', 'Damping',...
                        'XLim', XLimRidge, 'SquaredCWT', multiSignalMode);
                end
                if checkboxAmplDamp.Value % plot de l'amortissement
                    RidgeQtyPlot2(ridge, 'val', dampingRidgeName, 'EvaluationFunctionX', 'abs',...
                        'ScaleX', get(xscaleAmplDamp, 'String'), 'ScaleY', get(yscaleAmplDamp, 'String'),...
                        'Axes', axesFiguresCheckboxs2(5, kPlot), 'RenameAxes', true, 'NameX', 'Amplitude (m/s^2)', 'NameY', 'Damping',...
						'SquaredCWT', multiSignalMode);
                end
                if checkboxTimePhase.Value % plot de la phase
                    RidgeQtyPlot2(ridge, 'time', phaseRidgeName,...
                        'ScaleX', get(xscaleTimePhase, 'String'), 'ScaleY', get(yscaleTimePhase, 'String'),...
                        'Axes', axesFiguresCheckboxs2(6, kPlot), 'RenameAxes', true, 'NameX', 'Time (s)', 'NameY', 'Phase (rad)',...
                        'XLim', XLimRidge, 'SquaredCWT', multiSignalMode);
                end
            end
        end
        
        %% calcul des deformees
        if any([Checkboxs3.Value]) || any([Checkboxs4.Value])
            if multiSignalMode || autocorrelationMode
                if multiSignalMode
                    [timeShapes, freqsShapes, shapesShapes, amplitudesShapes] = ...
                        getModesSingleRidge(x(1,:), y, Q, fmin, fmax, NbFreq,...
                        'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu^2,...
                        'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects, 'MaxSlopeRidge', slopeRidge,...
                        'MotherWavelet', MotherWavelet);
                elseif autocorrelationMode
                    [timeShapesSVD, freqsShapesSVD, shapesShapesSVD, amplitudesShapesSVD] = ...
                        getModesCrossCorr(tRy, SVry, SVvectry, Q, fmin, fmax, NbFreq, autocorrelationNsvd,...
                        'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu,...
                        'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects, 'MaxSlopeRidge', slopeRidge,...
                        'MotherWavelet', MotherWavelet);
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
                        XquantLabel = 'Time (s)';
                        XquantScale = 'lin';
                        Xquantity = timeShapes;
                    elseif isequal(XquantName, 'ampl')
                        XquantLabel = 'Amplitude (m/s^2)';
                        XquantScale = 'lin';
                        Xquantity = absAmplitudesShapes;
                    elseif isequal(XquantName, 'log(ampl)')
                        XquantLabel = 'Amplitude (m/s^2)';
                        XquantScale = 'log';
                        Xquantity = absAmplitudesShapes;
                    elseif isequal(XquantName, 'freq')
                        XquantLabel = 'Frequence (Hz)';
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
                    end
                end
            else
                warning('! multi signal or autocorrelation modes only !');
            end
        end
        
        % changement du curseur
        set(fig, 'pointer', 'arrow');
        drawnow;
        
    end

%%

    function plotCrossCorr(Dx, y, NmaxLagCorr, forcePlot)
        % plot iif Ry change | forcePlot
        if nargin < 4
            forcePlot = isempty(tRy);
        end
        
        if isempty(tRy)
            tRy = Dx * (0:NmaxLagCorr);
            Ry = crossCorrelation(y, NmaxLagCorr);
        end
        
        if forcePlot
            figRy = figure('Name', ['crossCorr (', plotAxesName, ')']);
            axRy = axes(figRy);
            hold (axRy, 'on');
            
            Ndof = size(y, 1);
            legendRij = cell(1, Ndof^2);
            
            for i = 1:Ndof
                for j = 1:Ndof
                    if i == j
                        plot(axRy, tRy, reshape(Ry(i, j, :), [1, size(Ry, 3)]));
                    else
                        plot(axRy, tRy, reshape(Ry(i, j, :), [1, size(Ry, 3)]), ':');
                    end
                    legendRij{(i-1)*Ndof + j} = ['R', num2str(i), num2str(j)];
                end
            end
            
            xlim(axRy, [tRy(1), tRy(end)]);
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
        if isequal(hilbertPlotAxes, 0)
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

    function fourierTransform()
        getXY();
        fmin = eval(get(editfmin, 'String'));
        fmax = eval(get(editfmax, 'String'));
        
        if ~autocorrelationMode
            ffourier = figure;
            fourierPlotAxes = [];
            
            nbAxes = nbPlots;
            if autocorrelationMode
                nbAxes = nbPlots^2; % maximum
            end
            
            if multipleAxesDisplay %création des axes où sont plot les courbes
                for kPlot = 1:nbAxes
                    fourierPlotAxes(kPlot) = subplot0(nbAxes, 1, kPlot, axes(ffourier));
                    set(fourierPlotAxes(kPlot), 'XScale', 'lin', 'YScale', 'lin');
                end
            else
                axesffourier = axes(ffourier);
                hold(axesffourier, 'on');
                set(axesffourier, 'XScale', 'lin', 'YScale', 'lin');
                for kPlot = 1:nbAxes
                    fourierPlotAxes(kPlot) = axesffourier;
                end
            end
            
            if ~multiSignalMode % plot des courbes
                for kPlot = 1:nbAxes
                    Xfour = x(kPlot,:);
                    Yfour = y(kPlot,:);
                    Yfour = [Yfour, zeros(1, ZeroPaddingFourier*length(Yfour))];
                    Tfour = mean(diff(Xfour))*length(Yfour);
                    
                    four = fft(Yfour) / length(Yfour);
                    four = four(1:floor(end/2));
                    four(2:end) = 2*four(2:end);
                    freqs = 1/Tfour * (0:length(four)-1);
                    
                    hold(fourierPlotAxes(kPlot), 'on');
                    if isequal(FourierScale, 'lin')
                        plot(fourierPlotAxes(kPlot), freqs, abs(four));
                    elseif isequal(FourierScale, 'squared')
                        plot(fourierPlotAxes(kPlot), freqs, abs(four).^2);
                    elseif isequal(FourierScale, 'log')
                        plot(fourierPlotAxes(kPlot), freqs, abs(four));
                        set(fourierPlotAxes(kPlot), 'YScale', 'log');
                    elseif isequal(FourierScale, 'phase')
                        plot(fourierPlotAxes(kPlot), freqs, angle(four));
                    end
                    hold(fourierPlotAxes(kPlot), 'off');
                    
                    xlabel(fourierPlotAxes(kPlot), 'freq');
                    ylabel(fourierPlotAxes(kPlot), 'fft');
                    
                    %                 set(fourierPlotAxes(kPlot), 'Xlim', [fmin fmax]);
                end
            else
                FourierTot = 0;
                for kPlot = 1:nbAxes
                    Xfour = x(kPlot,:);
                    Yfour = y(kPlot,:);
                    Yfour = [Yfour, zeros(1, ZeroPaddingFourier*length(Yfour))];
                    Tfour = mean(diff(Xfour)) * length(Yfour);
                    four = fft(Yfour) / length(Yfour);
                    four = four(1:floor(end/2));
                    four(2:end) = 2*four(2:end);
                    FourierTot = FourierTot + four.^2;
                end
                freqs = 1/Tfour * (0:length(four)-1);
                four = sqrt(FourierTot);
                hold(fourierPlotAxes(kPlot), 'on');
                if isequal(FourierScale, 'lin')
                    plot(fourierPlotAxes(kPlot), freqs, abs(four));
                elseif isequal(FourierScale, 'squared')
                    plot(fourierPlotAxes(kPlot), freqs, abs(four).^2);
                elseif isequal(FourierScale, 'log')
                    plot(fourierPlotAxes(kPlot), freqs, abs(four));
                    set(fourierPlotAxes(kPlot), 'YScale', 'log');
                elseif isequal(FourierScale, 'phase')
                    plot(fourierPlotAxes(kPlot), freqs, angle(four));
                end
                hold(fourierPlotAxes(kPlot), 'off');
                xlabel(fourierPlotAxes(kPlot), 'freq');
                ylabel(fourierPlotAxes(kPlot), 'fft');
            end
            
        else % autocorr
            Dx = (x(1, end) - x(1, 1))/(size(x, 2)-1);
            NmaxLagCorr = floor(maxLagCorr/Dx);
            plotCrossCorr(Dx, y, NmaxLagCorr);
            
            if strcmp(autocorrelationFourierMenu.Checked, 'on')
                [SVfftrx, ~] = svdFFT(Ry, autocorrelationNsvd);
                
                for ksv = 1:autocorrelationNsvd
                    fourierFig = figure;
                    fourierPlotAxes = axes(fourierFig);
                    
                    Tfour = (tRy(end)-tRy(1)) / (length(tRy)-1) * length(tRy);
                    fftRx = SVfftrx{ksv}(1:floor(end/2));
                    freqs = 1/Tfour * (0:length(fftRx)-1);
                    
                    if isequal(FourierScale, 'lin')
                        plot(fourierPlotAxes, freqs, abs(fftRx));
                    elseif isequal(FourierScale, 'squared')
                        plot(fourierPlotAxes, freqs, abs(fftRx).^2);
                    elseif isequal(FourierScale, 'log')
                        plot(fourierPlotAxes, freqs, abs(fftRx));
                        set(fourierPlotAxes, 'YScale', 'log');
                    elseif isequal(FourierScale, 'phase')
                        plot(fourierPlotAxes, freqs, angle(fftRx));
                    end
                    xlabel(fourierPlotAxes, 'freq');
                    ylabel(fourierPlotAxes, 'fft');
                end
            else
                fourierFig = figure;
                fourierPlotAxes = axes(fourierFig);
                hold(fourierPlotAxes, 'on');
                
                Tfour = (tRy(end)-tRy(1)) / (length(tRy)-1) * length(tRy);
                
                legendRij = cell(1, nbPlots^2);
                for i = 1:nbPlots
                    for j = 1:nbPlots
                        fftRx = fft(reshape(Ry(i, j, :), [1, size(Ry, 3)]));
                        fftRx = fftRx(1:floor(end/2));
                        freqs = 1/Tfour * (0:length(fftRx)-1);
                        
                        if isequal(FourierScale, 'lin')
                            fourierPlot = plot(fourierPlotAxes, freqs, abs(fftRx));
                        elseif isequal(FourierScale, 'squared')
                            fourierPlot = plot(fourierPlotAxes, freqs, abs(fftRx).^2);
                        elseif isequal(FourierScale, 'log')
                            fourierPlot = plot(fourierPlotAxes, freqs, abs(fftRx));
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
                
                hold(fourierPlotAxes, 'off');
                xlabel(fourierPlotAxes, 'freq');
                ylabel(fourierPlotAxes, 'fft');
            end
        end
    end


buttonWavelet.Callback = @(~,~) show();

deleteButton.Callback = @(~,~) deletePlots();

buttonHilbert.Callback = @(~,~) hilbertTransform();

buttonFourier.Callback = @(~,~) fourierTransform();


end