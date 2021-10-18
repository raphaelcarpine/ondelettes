function [X, T] = removeNanSignal(X, T)
%REMOVENANSIGNAL Summary of this function goes here
%   Detailed explanation goes here

Nt0 = floor(size(X, 2)/2);
for Nt1 = Nt0:-1:1
    if any(isnan(X(:, Nt1)))
        Nt1 = Nt1 + 1;
        break
    end
end
for Nt2 = Nt0:size(X, 2)
    if any(isnan(X(:, Nt2)))
        Nt2 = Nt2 - 1;
        break
    end
end

if max(T(Nt1) - T(1), T(end) - T(Nt2)) > 0.1
    color = 2; % problem
else
    color = 1;
end

fprintf(color, 'Total sigal length: %.2fs ; NaN exclusion time : [%.2fs, %.2fs]\n',...
    [T(end)-T(1), T(Nt1) - T(1), T(Nt2) - T(end)]);

X = X(:, Nt1:Nt2);
T = T(Nt1:Nt2);

if ~isempty(T)
    T = T - T(1);
end

end

