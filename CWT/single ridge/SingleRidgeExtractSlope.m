function Fridge = SingleRidgeExtractSlope(T, freqs, CWT, slopeTimeConst)
%SINGLERIDGEEXTRACT Summary of this function goes here
%   Detailed explanation goes here

Nit = 2; % nb of iterations
Ndigitsmax = 5; % maximum of digits in the LongInt class base
precRidge = 1e-4; % precision of the ridge, i.e. ratio of the ridge impacted by the start and end values

%%

N = length(T);
dt = (T(end)-T(1))/(length(T)-1);
logFreqs = log(freqs);
M = abs(CWT).^2;

%% initial ridge extraction, lambda = 0

[M0, Phi] = max(M);
Phi = logFreqs(Phi);

Mref = mean(M0); % reference value for integral of M
lambda = Mref * slopeTimeConst^2;

prec = 2; % precision digits, in the LongInt base (2 corresponds to 64bit precision)
Prec = trunc(LongInt(1), -prec);

%% waitbar

[initWaitBar, updateWaitBar, closeWaitBar] = getWaitBar(N*Nit, 'displayTime', 0, ...
    'windowTitle', 'Computing ridge', 'msg', 'Initialization');
initWaitBar();

%% computation of ridge

for kit = 1:Nit
    updateWaitBar((kit-1)*N, sprintf('Iteration %d/%d', [kit, Nit]));
    
    [phiMax, Beta] = quadInterp(Phi); % computation of quadratic coefficients
    
    coeff = 2 + dt^2*Beta/lambda;
    coeffC = dt^2*Beta.*phiMax/lambda;
    
    Phi = nan(size(Phi));
    
    kt1 = 1; % subinterval computation start
    kt2 = 0; % previous subinterval end
    ktend = 1;
    while ktend < N
        % computation of ridge on subinterval
        
        % initialization
        A = LongInt.zeros(1, 2);
        B = LongInt.zeros(1, 2);
        C = LongInt.zeros(1, 2);
        A(1) = trunc(LongInt(1), -prec);
        B(2) = trunc(LongInt(1), -prec);
        
        % computation of phi_0 & phi_1 coefficients
        kt = 2;
        while kt1 + kt <= N
            coeffP = coeff(kt1+kt-1) .* Prec;
            A(kt+1) = trunc(coeffP .* A(kt), prec) - A(kt-1);
            B(kt+1) = trunc(coeffP .* B(kt), prec) - B(kt-1);
            C(kt+1) = trunc(coeffP .* C(kt), prec) - C(kt-1) - coeffC(kt1+kt-1) .* Prec;
            kt = kt + 1;
            if mod(kt, 10) == 0 && length(B(end).K) > Ndigitsmax
                break
            end
        end
        
        % computation of phi_0 & phi_{n-1} coefficients
        phiA = nan(size(A));
        phiB = nan(size(B));
        phiC = nan(size(C));
        for kt = 1:length(A)
            phiA(kt) = divideDouble(A(kt).*B(end) - B(kt).*A(end), B(end).trunc(-prec));
            phiB(kt) = divideDouble(B(kt), B(end));
            phiC(kt) = divideDouble(C(kt).*B(end) - B(kt).*C(end), B(end).trunc(-prec));
        end
        
        % computation of optimal phi_0 & phi_{n-1}
        dphiA = diff(phiA);
        dphiB = diff(phiB);
        dphiC = diff(phiC);
        beta = Beta(kt1:kt1+length(phiA)-1);
        HESS = [dt*sum(beta.*phiA.^2) + lambda/dt*sum(dphiA.^2), dt*sum(beta.*phiA.*phiB) + lambda/dt*sum(dphiA.*dphiB);
            dt*sum(beta.*phiA.*phiB) + lambda/dt*sum(dphiA.*dphiB), dt*sum(beta.*phiB.^2) + lambda/dt*sum(dphiB.^2)];
        R = [dt*sum(beta.*phiA.*(phiC-phiMax(kt1:kt1+length(phiA)-1))) + lambda/dt*sum(dphiA.*dphiC);
            dt*sum(beta.*phiB.*(phiC-phiMax(kt1:kt1+length(phiA)-1))) + lambda/dt*sum(dphiB.*dphiC)];
        phiEnds = - HESS\R;
        
        % computation of phi on subinterval, and integration to whole interval
        phiA = phiA(kt2-kt1+2:end);
        phiB = phiB(kt2-kt1+2:end);
        phiC = phiC(kt2-kt1+2:end);
        phi = phiA*phiEnds(1) + phiB*phiEnds(2) + phiC;
%         disp(phiA(1));
        Phi(kt2+1: kt2+length(phi)) = phi;
        ktend = kt2 + length(phi);
        kt1 = kt2 + sum(phiB < precRidge^2);
        kt2 = kt2 + sum(phiB < precRidge);
        
        % waitbar
        updateWaitBar((kit-1)*N + kt2);
    end
end

Fridge = exp(Phi);

closeWaitBar();


%% quadratic coefficients function

    function [phiMax, Beta] = quadInterp(logf)
        % closest frequencies indexes
        [~, Kf] = min(abs(logFreqs.' * ones(size(T)) - logf));
        Kf(Kf == 1) = 2;
        Kf(Kf == length(logFreqs)) = length(logFreqs) - 1;
        
        % M
        Mk = [M(Kf-1 + (0:size(M, 2)-1)*size(M, 1));
            M(Kf + (0:size(M, 2)-1)*size(M, 1));
            M(Kf+1 + (0:size(M, 2)-1)*size(M, 1))];
        logFk = [logFreqs(Kf-1); logFreqs(Kf); logFreqs(Kf+1)];
        
        % quad coeficients
        diffLogFk = [logFk(3, :)-logFk(2, :); logFk(1, :)-logFk(3, :); logFk(2, :)-logFk(1, :)];
        Beta = sum(Mk .* diffLogFk ./ prod(diffLogFk));
        phiMax = sum((Mk .* diffLogFk ./ prod(diffLogFk)).*(sum(logFk)-logFk) ./ (2*Beta));
    end

end

