function [X, T, DeltaT] = removeNanSignal(X, T)
%REMOVENANSIGNAL Summary of this function goes here
%   Detailed explanation goes here
%   DeltaT : décalage début signal

Nt0 = floor(size(X, 2)/4);
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
fprintf('NaN beginning: %.2fs\nNaN end: %.2fs\n', [T(Nt1) - T(1), T(end) - T(Nt2)]);
X = X(:, Nt1:Nt2);
T = T(:, Nt1:Nt2);
fprintf('Total sigal length: %.2fs\n', T(end)-T(1));

DeltaT = T(1);

T = T - T(1);

end

