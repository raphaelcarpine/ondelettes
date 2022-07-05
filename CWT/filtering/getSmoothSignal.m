function Xmean = getSmoothSignal(T, X, windowType, Tmean)
% windowType : 'rectangular', 'gaussian'
% Tmean : for 'rectangular' Tmean = b-a, for 'gaussian' Tmean = 2σ
dt = mean(diff(T));
if any(abs(diff(T)/dt - 1) > 1e-3) % pas de temps non constant
    error('uneven time step, not implemented');
else % pas de temps constant
    switch windowType
        case 'rectangular'
            halfNt = round(1/2*Tmean/dt);
            Nt = 2*halfNt + 1; % impair pour parité fenêtre
            localWindow = 1/Nt * ones(1, Nt);
        case 'gaussian'
            halfNt = round(2*Tmean/dt);
            Nt = 2*halfNt + 1; % impair pour parité fenêtre
            localWindow = exp(-(-halfNt:halfNt).^2 / (2*(Tmean/dt/2)^2));
            localWindow = localWindow / sum(localWindow);
    end
    
    Xmean = xcorr(X, localWindow);
    Xmean = Xmean(end-halfNt-length(X)+1 : end-halfNt);
    
    % correction bords
    localWindowSum = cumsum(localWindow);
    Xmean(1:halfNt+1) = Xmean(1:halfNt+1) ./ localWindowSum(halfNt+1:end);
    Xmean(end-halfNt:end) = Xmean(end-halfNt:end) ./ flip(localWindowSum(halfNt+1:end));
end
end