function ax = RidgeQtyPlot2(ridge, QtyX, QtyY ,varargin)
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

validQty = {'time', 'val', 'freq', 'diff', 'amor', 'freq2', 'pha'};
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

parse(p, ridge, QtyX, QtyY, varargin{:})

NameX = p.Results.NameX;
NameY = p.Results.NameY;
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

ax = p.Results.Axes;
if ax == 0
    fig = figure;
    ax = axes(fig);
end

%%
hold(ax, 'on');

for k_ridge = 1:length(ridge.freq)
    x = eval([EvaluationFunctionX, '(ridge.', QtyX, '{', num2str(k_ridge), '});']);
    y = eval([EvaluationFunctionY, '(ridge.', QtyY, '{', num2str(k_ridge), '});']);
    plot(x, y, 'Parent', ax);
    if showEdge
        %TODO
    end
end

%legend('show','location','best')
xlabel(ax ,NameX);
ylabel(ax, NameY);
ax.XScale = ScaleX;
ax.YScale = ScaleY;
ax.XGrid = 'on';
ax.YGrid = 'on';

hold(ax, 'off');
end