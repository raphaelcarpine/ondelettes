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
addParameter(p,'WaveletPlot', defaultWaveletPlot); %si les données viennent d'une courbe directement
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

if p.Results.WaveletPlot == 0
    getX = @() p.Results.X;
    getY = @() p.Results.Y;
else
    getX = @() get(p.Results.WaveletPlot, 'XData');
    getY = @() get(p.Results.WaveletPlot, 'YData');
end

fmin0 = p.Results.fmin;
fmax0 = p.Results.fmax;
NbFreq0 = p.Results.NbFreq;
Q0 = p.Results.Q;
MaxParallelRidges = p.Results.MaxParallelRidges;

%%
buttonWavelet = uicontrol('Parent',fig, 'Units', 'normalized','Style','pushbutton',...
    'String', 'ondelette');
paramPan = uipanel('Parent',fig, 'Units', 'normalized');
plotPan = uipanel('Parent',fig, 'Units', 'normalized');

buttonWavelet.Position = [0.02 0.02 0.96 0.15];
paramPan.Position = [0.02 0.19 0.4 0.79];
plotPan.Position = [0.44 0.19 0.54 0.79];

%%
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



%%
checkboxGeneral = uicontrol('Parent',plotPan, 'Units', 'normalized','Style','checkbox',...
    'String', 'general plot', 'Value', false);


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

Checkboxs = [checkboxGeneral, checkboxTimeAmpl, checkboxTimeFreq, checkboxAmplFreq];
XScales = [0, xscaleTimeAmpl, xscaleTimeFreq, xscaleAmplFreq];
YScales = [0, yscaleTimeAmpl, yscaleTimeFreq, yscaleAmplFreq];
n = length(Checkboxs);
for k=1:n
    Checkboxs(k).Position = [0.01, 0.01+(n-k)/n, 0.48, 1/n-0.02];
end
for k=2:n
    XScales(k).Position = [0.51, 0.01+(n-k)/n, 0.23, 1/n-0.02];
    YScales(k).Position = [0.76, 0.01+(n-k)/n, 0.23, 1/n-0.02];
end

scalesNames = {'linear', 'log'};
for k=2:n
    XScales(k).Callback = @(~,~) set(XScales(k), 'String', scalesNames{XScales(k).Value+1});
    YScales(k).Callback = @(~,~) set(YScales(k), 'String', scalesNames{YScales(k).Value+1});
end

LogScales = [yscaleTimeAmpl, xscaleAmplFreq];
for scale = LogScales
    set(scale, 'Value', 1);
    set(scale, 'String', scalesNames{scale.Value+1});
end



%%
    function show()
        fmin = str2double(get(editfmin, 'String'));
        fmax = str2double(get(editfmax, 'String'));
        NbFreq = str2double(get(editNbFreq, 'String'));
        Q = str2double(get(editQ, 'String'));
        PR = str2double(get(editPR, 'String')); %nombre max de ridges parallèles
        x = getX();
        y = getY();
        
        if checkboxGeneral.Value
            WvltPlot(x, y, linspace(fmin,fmax,NbFreq), Q);
        end
        
        ridge = RidgeExtract(x, y, Q, fmin, fmax, NbFreq, 'NbMaxParallelRidges', PR);
        
        if checkboxTimeAmpl.Value
            RidgeQtyPlot2(ridge, 'time', 'val', 'EvaluationFunctionY', 'abs',...
                'ScaleX', get(xscaleTimeAmpl, 'String'), 'ScaleY', get(yscaleTimeAmpl, 'String') );
        end
        if checkboxTimeFreq.Value
            RidgeQtyPlot2(ridge, 'time', 'freq',...
                'ScaleX', get(xscaleTimeFreq, 'String'), 'ScaleY', get(yscaleTimeFreq, 'String') );
        end
        if checkboxAmplFreq.Value
            RidgeQtyPlot2(ridge, 'val', 'freq', 'EvaluationFunctionX', 'abs',...
                'ScaleX', get(xscaleAmplFreq, 'String'), 'ScaleY', get(yscaleAmplFreq, 'String') );
        end
    end


buttonWavelet.Callback = @(~,~) show();


end