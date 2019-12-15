function fig = WaveletMenu(varargin)
%WaveletMenu Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

defaultX = nan;
defaultY = nan;
defaultFmin = 1;
defaultFmax = 10;
defaultNbFreq = 100;
defaultParent = 0;
defaultWaveletPlot = 0;
defaultQ = 1;
defaultMaxRidges = 1;
defaultMaxParallelRidges = 1;
defaultMultipleAxesDisplay = false;
defaultRidgeMinModu = 0;
defaultCtEdgeEffects = 3;
defaultZeroPaddingFourier = 0;
defaultMultiSignalMode = false;



checkParent = @(f) isa(f, 'matlab.ui.Figure') || isa(f, 'matlab.ui.container.Panel')...
    || isa(f, 'matlab.ui.container.Tab') || isa(f, 'matlab.ui.container.ButtonGroup');

addOptional(p, 'X', defaultX);
addOptional(p, 'Y', defaultY);
addParameter(p, 'fmin', defaultFmin);
addParameter(p, 'fmax', defaultFmax);
addParameter(p, 'NbFreq', defaultNbFreq);
addParameter(p,'Parent', defaultParent, checkParent);
addParameter(p,'WaveletPlot', defaultWaveletPlot); %si les données viennent d'une courbe directement (ou plusieurs)
addParameter(p,'Q', defaultQ);
addParameter(p,'MaxRidges', defaultMaxRidges);
addParameter(p,'MaxParallelRidges', defaultMaxParallelRidges);
addParameter(p,'MultipleAxesDisplay', defaultMultipleAxesDisplay);
addParameter(p,'RidgeMinModu', defaultRidgeMinModu);
addParameter(p,'CtEdgeEffects', defaultCtEdgeEffects);
addParameter(p,'ZeroPaddingFourier', defaultZeroPaddingFourier);
addParameter(p,'MultiSignalMode', defaultMultiSignalMode);

parse(p, varargin{:})

fig = p.Results.Parent;
if fig == 0
    fig = figure;
    fig.Units = 'characters';
    fig.Position(3) = 65;
    fig.Position(4) = 25;
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

if p.Results.WaveletPlot == 0
    getX = @() p.Results.X;
    getY = @() p.Results.Y;
    plotAxes = 0;
else
    getX = @() cellmat2mat (get(p.Results.WaveletPlot, 'XData'));
    getY = @() cellmat2mat (get(p.Results.WaveletPlot, 'YData'));
    plotAxes = cellmat2mat (get(p.Results.WaveletPlot, 'Parent'));
    
    for waveplt = p.Results.WaveletPlot
        parent = waveplt;
        while ~isa(parent, 'matlab.ui.Figure')
            parent = get(parent, 'Parent');
        end
        parent.CloseRequestFcn = @(~,~) closeParent(parent);
    end
end
nbPlots = size(getX(), 1);


fmin0 = p.Results.fmin;
fmax0 = p.Results.fmax;
NbFreq0 = p.Results.NbFreq;
Q0 = p.Results.Q;
MaxRidges = p.Results.MaxRidges;
MaxParallelRidges = p.Results.MaxParallelRidges;
multipleAxesDisplay = p.Results.MultipleAxesDisplay;
RidgeMinModu = p.Results.RidgeMinModu;
ctEdgeEffects = p.Results.CtEdgeEffects;
ZeroPaddingFourier = p.Results.ZeroPaddingFourier;
multiSignalMode = p.Results.MultiSignalMode;

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
paramPan = uipanel('Parent',waveletPan, 'Units', 'normalized');
plotPan = uipanel('Parent',waveletPan, 'Units', 'normalized');

buttonWavelet.Position = [0.02 0.02 0.96 0.13];
paramPan.Position = [0.02 0.16 0.4 0.81];
plotPan.Position = [0.44 0.16 0.54 0.81];

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
    'String', num2str(fmin0));

strfmax = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'fmax = ');
editfmax = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(fmax0));

strNbFreq = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'NbFreq = ');
editNbFreq = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(NbFreq0));

strQ = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'Q = ');
editQ = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(Q0));

strmaxR = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'max ridges');
editmaxR = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(MaxRidges));

strPR = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'max parallel ridges');
editPR = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(MaxParallelRidges));

Strs = [strfmin, strfmax, strNbFreq, strQ, strmaxR, strPR];
Edits =[editfmin, editfmax, editNbFreq, editQ, editmaxR, editPR];
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

checkboxTimeBand = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'time, bandwidth', 'Value', false);
xscaleTimeBand = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleTimeBand = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

checkboxAmplBand = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'ampl, bandwidth', 'Value', false);
xscaleAmplBand = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleAmplBand = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

checkboxTimePhase = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'time, phase', 'Value', false);
xscaleTimePhase = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);
yscaleTimePhase = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','togglebutton',...
    'String', 'linear', 'Value', false);

Checkboxs2 = [checkboxTimeAmpl, checkboxTimeFreq, checkboxAmplFreq, checkboxTimeBand, checkboxAmplBand,...
    checkboxTimePhase];
XScales = [xscaleTimeAmpl, xscaleTimeFreq, xscaleAmplFreq, xscaleTimeBand, xscaleAmplBand, xscaleTimePhase];
YScales = [yscaleTimeAmpl, yscaleTimeFreq, yscaleAmplFreq, yscaleTimeBand, yscaleAmplBand, yscaleTimePhase];
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




%% menu (autres paramètres)

%autres valeurs par défault
freqRidgeName = 'freq';
phaseRidgeName = 'pha2';

freqRidgeNames = {'freq', 'freq2'};
phaseRidgeNames = {'pha', 'pha2'};

fig.MenuBar = 'none';

paramMenu = uimenu(fig,'Text','Paramètres');

%freq
freqMenu = uimenu(paramMenu, 'Text','Frequence');
freqMenuChoices(1) = uimenu(freqMenu, 'Text', 'maximum module', 'Checked' ,'on');
freqMenuChoices(2) = uimenu(freqMenu, 'Text', 'derivée phase');
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
phaseMenuChoices(1) = uimenu(phaseMenu, 'Text', 'bornée');
phaseMenuChoices(2) = uimenu(phaseMenu, 'Text', 'continue', 'Checked' ,'on');
    function selectPhaseMenu(kchoice)
        for kchoices = 1:length(phaseMenuChoices)
            set(phaseMenuChoices(kchoices), 'Checked', 'off');
        end
        set(phaseMenuChoices(kchoice), 'Checked', 'on');
        
        phaseRidgeName = phaseRidgeNames{kchoice};
    end
set(phaseMenuChoices(1), 'CallBack', @(~,~) selectPhaseMenu(1));
set(phaseMenuChoices(2), 'CallBack', @(~,~) selectPhaseMenu(2));

%zero padding fourier
zeroPaddingFourierMenu = uimenu(paramMenu, 'Text','Zero padding Fourier');
    function setZeroPaddingFourier()
        ZeroPaddingFourier = nan;
        while isnan(ZeroPaddingFourier)
            ZeroPaddingFourier = inputdlg('Enter zero padding Fourier', 'zero padding Fourier');
            try
                ZeroPaddingFourier = str2double(ZeroPaddingFourier{1});
                ZeroPaddingFourier = round(ZeroPaddingFourier);
            catch
            end
        end
    end
set(zeroPaddingFourierMenu, 'CallBack', @(~,~) setZeroPaddingFourier);

%multipleAxesDisplay

multipleAxesDisplayMenu = uimenu(paramMenu, 'Text','multiple axes', 'Checked', multipleAxesDisplay);
    function switchMultipleAxesDisplay(~, ~)
        if strcmp(multipleAxesDisplayMenu.Checked, 'on')
            multipleAxesDisplayMenu.Checked = false;
            setMultipleAxesDisplay(false);
        else
            multipleAxesDisplayMenu.Checked = true;
            setMultipleAxesDisplay(true);
            multiSignalModeMenu.Checked = false;
            multiSignalMode = false;
        end
    end

multipleAxesDisplayMenu.MenuSelectedFcn = @switchMultipleAxesDisplay;

%multiSignalMode

multiSignalModeMenu = uimenu(paramMenu, 'Text','multi signal mode', 'Checked', multiSignalMode);
    function switchMultiSignalModeDisplay(~, ~)
        if strcmp(multiSignalModeMenu.Checked, 'on')
            multiSignalModeMenu.Checked = false;
            multiSignalMode = false;
        else
            multiSignalModeMenu.Checked = true;
            multipleAxesDisplayMenu.Checked = false;
            setMultipleAxesDisplay(false);
            multiSignalMode = true;
        end
    end

multiSignalModeMenu.MenuSelectedFcn = @switchMultiSignalModeDisplay;


%%

    function show()
        % evaluation des parametres
        fmin = eval(get(editfmin, 'String'));
        fmax = eval(get(editfmax, 'String'));
        NbFreq = eval(get(editNbFreq, 'String'));
        Q = eval(get(editQ, 'String'));
        maxR = eval(get(editmaxR, 'String')); %nombre max de ridges
        PR = eval(get(editPR, 'String')); %nombre max de ridges parallèles
        x = getX();
        y = getY();
        
        % plot de la transformee
        if checkboxGeneral.Value
            if ~multiSignalMode
                for kPlot = 1:nbPlots
                    WvltPlot(x(kPlot,:), y(kPlot,:), linspace(fmin,fmax,NbFreq), Q, 'ctEdgeEffects', ctEdgeEffects,...
                        'ctZeroPadding', ctEdgeEffects);
                end
            end
        end
        
        if checkboxModule.Value || checkboxPhase.Value
            if ~multiSignalMode
                for kPlot = 1:nbPlots
                    wavelet = WvltComp(x(kPlot,:), y(kPlot,:), linspace(fmin,fmax,NbFreq), Q, 'ct', ctEdgeEffects);
                    if checkboxModule.Value
                        WvltPlot2(x(kPlot,:), linspace(fmin,fmax,NbFreq), wavelet, 'module', Q, ctEdgeEffects);
                    end
                    if checkboxPhase.Value
                        WvltPlot2(x(kPlot,:), linspace(fmin,fmax,NbFreq), wavelet, 'phase', Q, ctEdgeEffects);
                    end
                end
            else
                wavelet = 0; % calcul de la somme des carrés de transformées
                for kPlot = 1:nbPlots
                    wavelet = wavelet +...
                        WvltComp(x(kPlot,:), y(kPlot,:), linspace(fmin,fmax,NbFreq), Q, 'ct', ctEdgeEffects).^2;
                end
                %wavelet = sqrt(wavelet);
                
                if checkboxModule.Value
                    WvltPlot2(x(kPlot,:), linspace(fmin,fmax,NbFreq), wavelet,...
                        'module', Q, ctEdgeEffects, 'sum wvlt^2');
                end
                if checkboxPhase.Value
                    WvltPlot2(x(kPlot,:), linspace(fmin,fmax,NbFreq), wavelet,...
                        'phase', Q, ctEdgeEffects, 'sum wvlt^2');
                end
            end
        end
        
        % calcul des ridges
        if ~multiSignalMode
            ridges = cell(1, nbPlots);
            for kPlot = 1:nbPlots
                ridges{kPlot} = RidgeExtract(x(kPlot,:), y(kPlot,:), Q, fmin, fmax, NbFreq,...
                    'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu,...
                    'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects);
            end
        else
            wavelet = 0;
            for kPlot = 1:nbPlots
                wavelet = wavelet +...
                    WvltComp(x(kPlot,:), y(kPlot,:), linspace(fmin,fmax,NbFreq), Q, 'ct', ctEdgeEffects).^2;
            end
            
            ridges = cell(1, 1);
            ridges{1} = RidgeExtract(x(kPlot,:), nan, Q, fmin, fmax, NbFreq,...
                'Wavelet', wavelet,...
                'NbMaxParallelRidges', PR, 'NbMaxRidges', maxR, 'MinModu', RidgeMinModu,...
                'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects);
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
        for kPlot = 1:length(ridges);
            ridge = ridges{kPlot};
            
            if ~isequal(plotAxes, 0) && checkboxTimeAmplPlot.Value % plot de l'amplitude directement sur l'axe
                newTimeAmplPlots = RidgeQtyPlot2(ridge, 'time', 'val', 'EvaluationFunctionY', 'abs',...
                    'Axes', plotAxes(kPlot), 'Grid', 'auto', 'RenameAxes', false);
                timeAmplPlots = [timeAmplPlots, newTimeAmplPlots];
            end
            
            if checkboxTimeAmpl.Value % plot de l'amplitude
                RidgeQtyPlot2(ridge, 'time', 'val', 'EvaluationFunctionY', 'abs',...
                    'ScaleX', get(xscaleTimeAmpl, 'String'), 'ScaleY', get(yscaleTimeAmpl, 'String'),...
                    'Axes', axesFiguresCheckboxs2(1, kPlot),...
                    'XLim', [x(kPlot,1), x(kPlot,end)]);
            end
            if checkboxTimeFreq.Value % plot de la frequence
                RidgeQtyPlot2(ridge, 'time', freqRidgeName,...
                    'ScaleX', get(xscaleTimeFreq, 'String'), 'ScaleY', get(yscaleTimeFreq, 'String'),...
                    'Axes', axesFiguresCheckboxs2(2, kPlot),...
                    'XLim', [x(kPlot,1), x(kPlot,end)]);
            end
            if checkboxAmplFreq.Value % plot de l'amplitude en fonction de la frequance
                RidgeQtyPlot2(ridge, 'val', freqRidgeName, 'EvaluationFunctionX', 'abs',...
                    'ScaleX', get(xscaleAmplFreq, 'String'), 'ScaleY', get(yscaleAmplFreq, 'String'),...
                    'Axes', axesFiguresCheckboxs2(3, kPlot));
            end
            if checkboxTimeBand.Value % plot de l'amortissement
                RidgeQtyPlot2(ridge, 'time', 'bandwidth',...
                    'ScaleX', get(xscaleTimeBand, 'String'), 'ScaleY', get(yscaleTimeBand, 'String'),...
                    'Axes', axesFiguresCheckboxs2(4, kPlot),...
                    'XLim', [x(kPlot,1), x(kPlot,end)]);
            end
            if checkboxAmplBand.Value % plot de l'amortissement
                RidgeQtyPlot2(ridge, 'val', 'bandwidth', 'EvaluationFunctionX', 'abs',...
                    'ScaleX', get(xscaleAmplBand, 'String'), 'ScaleY', get(yscaleAmplBand, 'String'),...
                    'Axes', axesFiguresCheckboxs2(5, kPlot));
            end
            if checkboxTimePhase.Value % plot de la phase
                RidgeQtyPlot2(ridge, 'time', phaseRidgeName,...
                    'ScaleX', get(xscaleTimePhase, 'String'), 'ScaleY', get(yscaleTimePhase, 'String'),...
                    'Axes', axesFiguresCheckboxs2(6, kPlot),...
                    'XLim', [x(kPlot,1), x(kPlot,end)]);
            end
        end
    end




    function deletePlots()
        try
            delete(timeAmplPlots);
        catch
        end
        timeAmplPlots = [];
    end





    function hilbertTransform()
        x = getX();
        y = getY();
        
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





    function fourierTransform()
        x = getX();
        y = getY();
        fmin = eval(get(editfmin, 'String'));
        fmax = eval(get(editfmax, 'String'));
        
        ffourier = figure;
        fourierPlotAxes = [];
        if multipleAxesDisplay %création des axes où sont plot les courbes
            for kPlot = 1:nbPlots
                fourierPlotAxes(kPlot) = subplot0(nbPlots, 1, kPlot, axes(ffourier));
                set(fourierPlotAxes(kPlot), 'XScale', 'lin', 'YScale', 'lin');
            end
        else
            axesffourier = axes(ffourier);
            hold(axesffourier, 'on');
            set(axesffourier, 'XScale', 'lin', 'YScale', 'lin');
            for kPlot = 1:nbPlots
                fourierPlotAxes(kPlot) = axesffourier;
            end            
        end
        
        if ~multiSignalMode % plot des courbes
            for kPlot = 1:nbPlots
                Xfour = x(kPlot,:);
                Yfour = y(kPlot,:);
                Yfour = [Yfour, zeros(1, ZeroPaddingFourier*length(Yfour))];
                Tfour = (Xfour(end)-Xfour(1))*length(Yfour)/length(Xfour);
                
                four = fft(Yfour);
                four = four(1:floor(end/2));
                freqs = linspace(0, length(four)/Tfour, length(four));
                
                hold(fourierPlotAxes(kPlot), 'on');
                plot(fourierPlotAxes(kPlot), freqs, abs(four));
                hold(fourierPlotAxes(kPlot), 'off');
                
                xlabel(fourierPlotAxes(kPlot), 'freq');
                ylabel(fourierPlotAxes(kPlot), 'fft');
%                 set(fourierPlotAxes(kPlot), 'Xlim', [fmin fmax]);
            end
        else
            FourierTot = 0;
            for kPlot = 1:nbPlots
                Xfour = x;
                Yfour = y(kPlot,:);
                Yfour = [Yfour, zeros(1, ZeroPaddingFourier*length(Yfour))];
                Tfour = (Xfour(end)-Xfour(1))*length(Yfour)/length(Xfour);
                four = fft(Yfour);
                four = four(1:floor(end/2));
                FourierTot = FourierTot + four.^2;
            end
            freqs = linspace(0, length(four)/Tfour, length(four));
            hold(fourierPlotAxes(kPlot), 'on');
            plot(freqs, sqrt(abs(FourierTot)), 'Parent', fourierPlotAxes(kPlot));
            hold(fourierPlotAxes(kPlot), 'off');
            xlabel(fourierPlotAxes(kPlot), 'freq');
            ylabel(fourierPlotAxes(kPlot), 'fft');
        end
        
    end


buttonWavelet.Callback = @(~,~) show();

deleteButton.Callback = @(~,~) deletePlots();

buttonHilbert.Callback = @(~,~) hilbertTransform();

buttonFourier.Callback = @(~,~) fourierTransform();


end