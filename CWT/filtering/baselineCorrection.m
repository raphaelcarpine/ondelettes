function Xbaseline = baselineCorrection(T, X, n)
%BASELINECORRECTION Summary of this function goes here
%   Detailed explanation goes here

Tpuiss = ones(n+1, length(T));
for k = 1:n
    Tpuiss(k+1, :) = T.^k;
end

c = Tpuiss.' \ X.';

Xbaseline = X - c.'*Tpuiss;

end

