function meanQty = RidgeQtyAverage(ridge, Qty, varargin)
%% Quantity : time, val, freq, diff, amor, freq2, pha
p = inputParser;

defaultWeightedByAmpl = false;
defaultEvaluationFunction = '';

validQty = {'val', 'freq', 'diff', 'damping', 'damping2', 'damping3', 'bandwidth', 'freq2', 'pha', 'pha2'};
checkQty = @(str) ismember(str, validQty);
validEvaluationFunction = {'', 'abs', 'angle', 'real', 'imag', 'log'};
checkEvaluationFunction = @(str) ismember(str, validEvaluationFunction);

addRequired(p, 'ridge');
addRequired(p, 'Qty', checkQty);
addParameter(p,'WeightedByAmpl', defaultWeightedByAmpl);
addParameter(p,'EvaluationFunction', defaultEvaluationFunction, checkEvaluationFunction);

parse(p, ridge, Qty, varargin{:})

WeightedByAmpl = p.Results.WeightedByAmpl;
EvaluationFunction = p.Results.EvaluationFunction;


%% 

meanQty = nan(size(ridge.freq));

for k_ridge = 1:length(ridge.freq)
    ampl = abs(ridge.val{k_ridge});
    qty = eval([EvaluationFunction, '(ridge.', Qty, '{', num2str(k_ridge), '});']);
    
    if WeightedByAmpl
        weights = ampl / sum(ampl);
    else
        weights = ones(size(ampl)) / length(ampl);
    end
    
    meanQty(k_ridge) = mean(qty.*weights);
end

end