function fig = WaveletMenu(fmin,fmax,NbFreq, varargin)
%WaveletMenu Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

defaultX = nan;
defaultY = nan;
defaultParent = 0;
defaultWaveletPlot = 0;

checkParent = @(f) isa(f, 'matlab.ui.Figure') || isa(f, 'matlab.ui.container.Panel')...
    || isa(f, 'matlab.ui.container.Tab') || isa(f, 'matlab.ui.container.ButtonGroup');

addRequired(p, 'fmin');
addRequired(p, 'fmax');
addRequired(p, 'NbFreq');
addOptional(p, 'X', defaultX);
addOptional(p, 'Y', defaultY);
addParameter(p,'Parent', defaultParent, checkParent);
addParameter(p,'WaveletPlot', defaultWaveletPlot);

parse(p, fmin,fmax,NbFreq, varargin{:})

fig = p.Results.Parent;
if fig == 0
    fig = figure;
    fig.Units = 'characters';
    fig.Position(3) = 34;
    fig.Position(4) = 11;
    fig.MenuBar = 'none';
%     fig.ToolBar = 'none';
end

if p.Results.WaveletPlot == 0
    getX = @() X;
    getY = @() Y;
else
    getX = @() get(p.Results.WaveletPlot, 'XData');
    getY = @() get(p.Results.WaveletPlot, 'YData');
end

%%
buttonWavelet = uicontrol('Parent',fig, 'Units', 'characters','Position', [21 0.5 12 1.5],'Style','pushbutton',...
    'String', 'ondelette');

textQ = uicontrol('Parent',fig, 'Units', 'characters','Position', [1 0.5 5 1.5],'Style','text',...
    'String', 'Q = ');
editQ = uicontrol('Parent',fig, 'Units', 'characters','Position', [6 0.5 12 1.5],'Style','edit',...
    'String', '1');

textPR = uicontrol('Parent',fig, 'Units', 'characters','Position', [1 2.5 23 1.5],'Style','text',...
    'String', 'max parallel ridges : ');
editPR = uicontrol('Parent',fig, 'Units', 'characters','Position', [23 2.5 8 1.5],'Style','edit',...
    'String', '1');

checkboxAmplFreq = uicontrol('Parent',fig, 'Units', 'characters','Position', [1 5 32 1],'Style','checkbox',...
    'String', 'amplitude, frequence', 'Value', false);
checkboxTimeFreq = uicontrol('Parent',fig, 'Units', 'characters','Position', [1 6.5 32 1],'Style','checkbox',...
    'String', 'time, frequence', 'Value', true);
checkboxTimeAmpl = uicontrol('Parent',fig, 'Units', 'characters','Position', [1 8 32 1],'Style','checkbox',...
    'String', 'time, amplitude', 'Value', true);
checkboxGeneral = uicontrol('Parent',fig, 'Units', 'characters','Position', [1 9.5 32 1],'Style','checkbox',...
    'String', 'general plot', 'Value', true);

%%
    function show()
        Q = str2double(get(editQ, 'String'));
        PR = str2double(get(editPR, 'String')); %nombre max de ridges parallèles
        x = getX();
        y = getY();
        
        if checkboxGeneral.Value
            WvltPlot(x, y, linspace(fmin,fmax,NbFreq), Q);
        end
        
        ridge = RidgeExtract(x, y, Q, fmin, fmax, NbFreq, 'NbMaxParallelRidges', PR);
        
        if checkboxTimeAmpl.Value
            RidgeQtyPlot2(ridge, 'time', 'val', 'EvaluationFunctionY', 'abs', 'ScaleY', 'log')
        end
        if checkboxTimeFreq.Value
            RidgeQtyPlot2(ridge, 'time', 'freq')
        end
        if checkboxAmplFreq.Value
            RidgeQtyPlot2(ridge, 'val', 'freq', 'EvaluationFunctionX', 'abs', 'ScaleX', 'log')
        end
    end


buttonWavelet.Callback = @(~,~) show();


end