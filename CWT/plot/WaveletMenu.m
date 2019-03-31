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
defaultMaxParallelRidges = 1;

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
addParameter(p,'MaxParallelRidges', defaultMaxParallelRidges);

parse(p, varargin{:})

fig = p.Results.Parent;
if fig == 0
    fig = figure;
    fig.Units = 'characters';
    fig.Position(3) = 65;
    fig.Position(4) = 15;
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
        parent.CloseRequestFcn = @(~,~) delete([parent, fig]);
    end
end
nbPlots = size(getX(), 1);


fmin0 = p.Results.fmin;
fmax0 = p.Results.fmax;
NbFreq0 = p.Results.NbFreq;
Q0 = p.Results.Q;
MaxParallelRidges = p.Results.MaxParallelRidges;


%% bouton ondelettes et panneaux param et sorties
buttonWavelet = uicontrol('Parent',fig, 'Units', 'normalized','Style','pushbutton',...
    'String', 'wavelet');
paramPan = uipanel('Parent',fig, 'Units', 'normalized');
plotPan = uipanel('Parent',fig, 'Units', 'normalized');

buttonWavelet.Position = [0.02 0.02 0.96 0.15];
paramPan.Position = [0.02 0.19 0.4 0.79];
plotPan.Position = [0.44 0.19 0.54 0.79];

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

strPR = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','text',...
    'String', 'max parallel ridges : ');
editPR = uicontrol('Parent',paramPan, 'Units', 'normalized','Style','edit',...
    'String', num2str(MaxParallelRidges));

Strs = [strfmin, strfmax, strNbFreq, strQ, strPR];
Edits =[editfmin, editfmax, editNbFreq, editQ, editPR];
n = length(Strs);
for k=1:n
    Strs(k).Position = [0.01, 0.01+(n-k)/n, 0.48, 1/n-0.02];
    Edits(k).Position = [0.51, 0.01+(n-k)/n, 0.48, 1/n-0.02];
end



%% sorties

checkboxGeneral = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'general plot', 'Value', false);

Checkboxs1 = [checkboxGeneral];

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

Checkboxs2 = [checkboxTimeAmpl, checkboxTimeFreq, checkboxAmplFreq];
XScales = [xscaleTimeAmpl, xscaleTimeFreq, xscaleAmplFreq];
YScales = [yscaleTimeAmpl, yscaleTimeFreq, yscaleAmplFreq];
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



%%
    function show()
        fmin = eval(get(editfmin, 'String'));
        fmax = eval(get(editfmax, 'String'));
        NbFreq = eval(get(editNbFreq, 'String'));
        Q = eval(get(editQ, 'String'));
        PR = eval(get(editPR, 'String')); %nombre max de ridges parallèles
        x = getX();
        y = getY();
        
        if checkboxGeneral.Value
            for kPlot = 1:nbPlots
                WvltPlot(x(kPlot,:), y(kPlot,:), linspace(fmin,fmax,NbFreq), Q);
            end
        end
        
        ridges = {};
        for kPlot = 1:nbPlots
            ridges{end+1} = RidgeExtract(x(kPlot,:), y(kPlot,:), Q, fmin, fmax, NbFreq, 'NbMaxParallelRidges', PR);
        end
        
        for kCheck = 1:length(Checkboxs2)
            cb = Checkboxs2(kCheck);
            if cb.Value
                FiguresCheckboxs2(kCheck) = figure;
            end
        end
        for kPlot = 1:nbPlots
            ridge = ridges{kPlot};
            
            if ~isequal(plotAxes, 0) && checkboxTimeAmplPlot.Value
                newTimeAmplPlots = RidgeQtyPlot2(ridge, 'time', 'val', 'EvaluationFunctionY', 'abs',...
                    'Axes', plotAxes(kPlot), 'Grid', 'auto');
                timeAmplPlots = [timeAmplPlots, newTimeAmplPlots];
            end
            
            if checkboxTimeAmpl.Value
                RidgeQtyPlot2(ridge, 'time', 'val', 'EvaluationFunctionY', 'abs',...
                    'ScaleX', get(xscaleTimeAmpl, 'String'), 'ScaleY', get(yscaleTimeAmpl, 'String'),...
                    'Axes', subplot(nbPlots, 1, kPlot, axes(FiguresCheckboxs2(1))),...
                    'XLim', [x(kPlot,1), x(kPlot,end)]);
            end
            if checkboxTimeFreq.Value
                RidgeQtyPlot2(ridge, 'time', 'freq',...
                    'ScaleX', get(xscaleTimeFreq, 'String'), 'ScaleY', get(yscaleTimeFreq, 'String'),...
                    'Axes', subplot(nbPlots, 1, kPlot, axes(FiguresCheckboxs2(2))),...
                    'XLim', [x(kPlot,1), x(kPlot,end)]);
            end
            if checkboxAmplFreq.Value
                RidgeQtyPlot2(ridge, 'val', 'freq', 'EvaluationFunctionX', 'abs',...
                    'ScaleX', get(xscaleAmplFreq, 'String'), 'ScaleY', get(yscaleAmplFreq, 'String'),...
                    'Axes', subplot(nbPlots, 1, kPlot, axes(FiguresCheckboxs2(3))));
            end
        end
    end

    function deletePlots()
        delete(timeAmplPlots);
        timeAmplPlots = [];
    end


buttonWavelet.Callback = @(~,~) show();

deleteButton.Callback = @(~,~) deletePlots();


end