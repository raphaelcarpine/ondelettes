function [A0, coeff] = regLinFunc(data)
%REGLINFUNC Summary of this function goes here
%   Detailed explanation goes here

useArea = false;

if useArea
    for k0 = 1:length(data)
        [~, I] = sort(data{k0}(1, :));
        data{k0} = data{k0}(:, I);
    end
    
    S = @(p) sum(cellfun(@(d) sum(diff(d(1, :)) .*...
        ((d(2, 1:end-1) - (p(1)+p(2)*d(1, 1:end-1))).^2 + (d(2, 2:end) - (p(1)+p(2)*d(1, 2:end))).^2)/2),...
        data));
    
    optionsReg = optimoptions(@lsqnonlin, 'MaxIterations', 1e8,...
        'StepTolerance', 1e-10, 'MaxFunctionEvaluations', inf, 'FunctionTolerance', 0);
    
    P = lsqnonlin(S, [0 0], [], [], optionsReg);
else
    data = cell2mat(data);
    X = data(1, :).';
    Y = data(2, :).';
    P = [ones(size(X)), X] \ Y;
end

A0 = P(1);
coeff = P(2);

end

