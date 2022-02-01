function Xmean = getSmoothSignal(T, X, windowType, Tmean)
dt = mean(diff(T));
if any(abs(diff(T)/dt - 1) > 1e-4) % pas de temps non constant
    error('uneven time step, not implemented');
else % pas de temps constant
    halfNt = round(1/2*Tmean/dt);
    Nt = 2*halfNt + 1; % impair pour parité fenêtre
    switch windowType
        case 'rectangular'
            localWindow = 1/Nt * ones(1, Nt);
        case 'gaussian'
            localWindow = exp(-(-halfNt:halfNt).^2 / (2*(halfNt/3)^2));
            localWindow = localWindow / sum(localWindow);
    end
    
    Xmean = xcorr(X, localWindow);
    Xmean = Xmean(end-halfNt-length(X)+1 : end-halfNt);
end
end