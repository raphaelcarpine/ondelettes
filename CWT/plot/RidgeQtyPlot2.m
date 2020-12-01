function plts = RidgeQtyPlot2(ridge, QtyX, QtyY ,varargin)
%% Quantity : time, val, freq, diff, amor, freq2, pha
p = inputParser;

defaultNameX = QtyX;
defaultNameY = QtyY;
defaultAxes = 0;
defaultScaleX = 'linear';
defaultScaleY = 'linear';
defaultEvaluationFunctionX = '';
defaultEvaluationFunctionY = '';
defaultShowEdge = true;
defaultXLim = nan;
defaultYLim = nan;
defaultThreshold = 0;
defaultGrid = 'on'; % 'auto' pour ne pas changer
defaultRenameAxes = true;

validQty = {'time', 'val', 'freq', 'diff', 'damping', 'damping2', 'damping3', 'bandwidth', 'freq2', 'pha', 'pha2'};
checkQty = @(str) ismember(str, validQty);
validEvaluationFunction = {'', 'abs', 'angle', 'real', 'imag', 'log'};
checkEvaluationFunction = @(str) ismember(str, validEvaluationFunction);

addRequired(p, 'ridge');
addRequired(p, 'QtyX', checkQty);
addRequired(p, 'QtyY', checkQty);
addParameter(p, 'NameX', defaultNameX, @ischar);
addParameter(p, 'NameY', defaultNameY, @ischar);
addParameter(p,'Axes', defaultAxes);
addParameter(p,'ScaleX', defaultScaleX);
addParameter(p,'ScaleY', defaultScaleY);
addParameter(p,'EvaluationFunctionX', defaultEvaluationFunctionX, checkEvaluationFunction);
addParameter(p,'EvaluationFunctionY', defaultEvaluationFunctionY, checkEvaluationFunction);
addParameter(p,'showEdge', defaultShowEdge);
addParameter(p,'XLim', defaultXLim);
addParameter(p,'YLim', defaultYLim);
addParameter(p,'Threshold', defaultThreshold);
addParameter(p,'Grid', defaultGrid);
addParameter(p,'RenameAxes', defaultRenameAxes);

parse(p, ridge, QtyX, QtyY, varargin{:})

NameX = p.Results.NameX;
NameY = p.Results.NameY;
if isequal(NameX, 'val')
    NameX = 'CWT';
end
if isequal(NameY, 'val')
    NameY = 'CWT';
end
EvaluationFunctionX = p.Results.EvaluationFunctionX;
EvaluationFunctionY = p.Results.EvaluationFunctionY;
% if ~isequal(EvaluationFunctionX, '')
%     NameX = [EvaluationFunctionX, '(', NameX, ')'];
% end
% if ~isequal(EvaluationFunctionY, '')
%     NameY = [EvaluationFunctionY, '(', NameY, ')'];
% end
ScaleX = p.Results.ScaleX;
ScaleY = p.Results.ScaleY;

showEdge = p.Results.showEdge;
XLim = p.Results.XLim;
YLim = p.Results.YLim;
Threshold = p.Results.Threshold;
Grid = p.Results.Grid;

ax = p.Results.Axes;
if isequal(ax, 0)
    fig = figure;
    ax = axes(fig);
end

RenameAxes = p.Results.RenameAxes;

%% 

hold(ax, 'on');

plts = [];
for k_ridge = 1:length(ridge.freq)
    x = eval([EvaluationFunctionX, '(ridge.', QtyX, '{', num2str(k_ridge), '});']);
    y = eval([EvaluationFunctionY, '(ridge.', QtyY, '{', num2str(k_ridge), '});']);
    
    plts(end+1) = plot(x, y, 'Parent', ax);
    if showEdge
        %TODO
    end
end

if strcmp(QtyY, 'val') && Threshold > 0
    yline(ax, Threshold, 'r--');
    thresholdTxt = text(ax, ax.XLim(1), Threshold, ' threshold', 'VerticalAlignment', 'bottom', 'Color', 'r');
    addlistener(ax, 'XLim', 'PostSet', @(~,~) set(thresholdTxt, 'Position', [ax.XLim(1), Threshold]));
end
    

%legend('show','location','best')
if RenameAxes
    xlabel(ax ,NameX);
    ylabel(ax, NameY);
end
ax.XScale = ScaleX;
ax.YScale = ScaleY;

if ~isequal(Grid, 'auto')
    ax.XGrid = Grid;
    ax.YGrid = Grid;
end

if ~all(isnan(XLim))
    ax.XLim = XLim;
end
if ~all(isnan(YLim))
    ax.YLim = YLim;
end

hold(ax, 'off');
end