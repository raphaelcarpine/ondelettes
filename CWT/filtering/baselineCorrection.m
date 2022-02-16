function Xbaseline = baselineCorrection(T, X, n)
%BASELINECORRECTION Summary of this function goes here
%   Detailed explanation goes here

Tnorm = T - T(1);
Tnorm = Tnorm / Tnorm(end);

Tpuiss = ones(n+1, length(Tnorm));
for k = 1:n
    Tpuiss(k+1, :) = Tnorm.^k;
end

c = Tpuiss.' \ X.';

Xbaseline = X - c.'*Tpuiss;

end

