function [trdt, Xrdt] = getRDT(t, X, Tmax, KTrdt)
%GETRDT Random Decrement Technique
%   


%%

dt = mean(diff(t(1, :)));
NTmax = floor(Tmax / dt) + 1;
if NTmax == 1
    warning('Tmax < dt');
end


%% RDT

trdt = dt * (0:NTmax-1);

Xrdt = zeros(size(X, 1), NTmax);
for kTrdt = KTrdt
    Xrdt = Xrdt + X(:, kTrdt:kTrdt+NTmax-1);
end
Xrdt = Xrdt / length(KTrdt);

end

