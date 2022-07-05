function [X, T] = removeRedundantData(X, T)
% enlever donnees redondantes (double t à cause du mode transmit/log)
%   Detailed explanation goes here

dt0 = mean(diff(T));
T0 = T;
X0 = X;
T = nan(size(T));
X = nan(size(X));

kt = 1;
kt0 = 1;
while kt0 <= length(T0)-1
    if abs(T0(kt0) - T0(kt0+1)) < 1e-2 * dt0
        if all(~isnan(X0(:, kt0)))
            X(:, kt) = X0(:, kt0);
        elseif all(~isnan(X0(:, kt0+1)))
            X(:, kt) = X0(:, kt0+1);
        else
            for kc = 1:size(X0, 1)
                if ~isnan(X0(kc, kt0))
                    X(kc, kt) = X0(kc, kt0);
                else
                    X(kc, kt) = X0(kc, kt0+1);
                end
            end
%             if any(~isnan(X0(:, kt0:kt0+1)), 'all') % debug
%                 disp(X0(:, kt0:kt0+1));
%             end
%             if any(~isnan(X0(:, kt0:kt0+1)), 'all') && any(isnan(X0(:, kt0)) & isnan(X0(:, kt0+1))) % debug
%                 disp(X0(:, kt0:kt0+1));
%             end
        end
        T(kt) = T0(kt0);
        kt0 = kt0+2;
    else
        X(:, kt) = X0(:, kt0);
        T(kt) = T0(kt0);
        kt0 = kt0+1;
    end
    kt = kt+1;
end
T = T(1:kt-1);
X = X(:, 1:kt-1);

%% test

if any(abs(diff(T)/mean(diff(T))-1) > 1e-2)
    error('time error');
end

end

