function Fridge = SingleRidgeExtractSlope(T, freqs, CWT, slopeTimeConst)
%SINGLERIDGEEXTRACT Summary of this function goes here
%   Detailed explanation goes here

N = length(T);
dt = (T(end)-T(1))/(length(T)-1);
logFreqs = log(freqs);
M = abs(CWT).^2;

%% initial ridge extraction, lambda = 0

[M0, phi0] = max(M);
phi0 = logFreqs(phi0);

Mref = mean(M0); % reference value for integral of M
lambda = Mref * slopeTimeConst^2;

%% computation of quadratic coefficients

[phi0, Beta] = quadInterp(phi0);

coeff = 2 + dt^2*Beta/lambda;
coeffC = dt^2*Beta.*phi0/lambda;

prec = 2; % precision digits, in the LongInt base (2 corresponds to a 64bit precision)
Prec = trunc(LongInt(1), -prec);

A = LongInt.zeros(1, N);
B = LongInt.zeros(1, N);
C = LongInt.zeros(1, N);
A(1) = trunc(LongInt(1), -prec);
B(2) = trunc(LongInt(1), -prec);

for kt = 3:N
    A(kt) = trunc(coeff(kt-1) .* Prec .* A(kt-1), prec) - A(kt-2);
    B(kt) = trunc(coeff(kt-1) .* Prec .* B(kt-1), prec) - B(kt-2);
    C(kt) = trunc(coeff(kt-1) .* Prec .* C(kt-1), prec) - C(kt-2) - coeffC(kt-1) .* Prec;
end

phiA = nan(1, N);
phiB = nan(1, N);
phiC = nan(1, N);
for kt = 1:N
    phiA(kt) = divideDouble(A(kt).*B(end) - B(kt).*A(end), B(end).trunc(-prec));
    phiB(kt) = divideDouble(B(kt), B(end));
    phiC(kt) = divideDouble(C(kt).*B(end) - B(kt).*C(end), B(end).trunc(-prec));
end

dphiA = diff(phiA);
dphiB = diff(phiB);
dphiC = diff(phiC);


HESS = [dt*sum(Beta.*phiA.^2) + lambda/dt*sum(dphiA.^2), dt*sum(Beta.*phiA.*phiB) + lambda/dt*sum(dphiA.*dphiB);
    dt*sum(Beta.*phiA.*phiB) + lambda/dt*sum(dphiA.*dphiB), dt*sum(Beta.*phiB.^2) + lambda/dt*sum(dphiB.^2)];
R = [dt*sum(Beta.*phiA.*(phiC-phi0)) + lambda/dt*sum(dphiA.*dphiC);
    dt*sum(Beta.*phiB.*(phiC-phi0)) + lambda/dt*sum(dphiB.*dphiC)];
phiEnds = - HESS\R;

phi = phiA*phiEnds(1) + phiB*phiEnds(2) + phiC;

Fridge = exp(phi);

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

