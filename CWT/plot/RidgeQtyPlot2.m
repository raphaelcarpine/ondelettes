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
defaultGrid = 'on'; % 'auto' pour ne pas changer
defaultRenameAxes = true;
defaultSquaredCWT = false;

validQty = {'time', 'val', 'freq', 'diff', 'damping', 'bandwidth', 'freq2', 'pha', 'pha2'};
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
addParameter(p,'Grid', defaultGrid);
addParameter(p,'RenameAxes', defaultRenameAxes);
addParameter(p,'SquaredCWT', defaultSquaredCWT);

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
if ~isequal(EvaluationFunctionX, '')
    NameX = [EvaluationFunctionX, '(', NameX, ')'];
end
if ~isequal(EvaluationFunctionY, '')
    NameY = [EvaluationFunctionY, '(', NameY, ')'];
end
ScaleX = p.Results.ScaleX;
ScaleY = p.Results.ScaleY;

showEdge = p.Results.showEdge;
XLim = p.Results.XLim;
YLim = p.Results.YLim;
Grid = p.Results.Grid;

ax = p.Results.Axes;
if isequal(ax, 0)
    fig = figure;
    ax = axes(fig);
end

RenameAxes = p.Results.RenameAxes;
SquaredCWT = p.Results.SquaredCWT;

%% 

hold(ax, 'on');

plts = [];
for k_ridge = 1:length(ridge.freq)
    x = eval([EvaluationFunctionX, '(ridge.', QtyX, '{', num2str(k_ridge), '});']);
    y = eval([EvaluationFunctionY, '(ridge.', QtyY, '{', num2str(k_ridge), '});']);
    
    if SquaredCWT
        if isequal(QtyX, 'val')
            x = sqrt(x);
        elseif ismember(QtyX, {'damping', 'bandwidth', 'freq2', 'pha', 'pha2'})
            x = x/2;
        end
        if isequal(QtyY, 'val')
            y = sqrt(y);
        elseif ismember(QtyY, {'damping', 'bandwidth', 'freq2', 'pha', 'pha2'})
            y = y/2;
        end
    end
    
    plts(end+1) = plot(x, y, 'Parent', ax);
    if showEdge
        %TODO
    end
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