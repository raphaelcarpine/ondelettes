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

[initWaitBar, updateWaitBar, closeWaitBar] =...
    getWaitBar(length(KTrdt), 'displayTime', 5, 'windowTitle', 'Computing random decrement');
initWaitBar();

trdt = dt * (0:NTmax-1);
Xrdt = zeros(size(X, 1), NTmax);
for i_k = 1:length(KTrdt)
    kTrdt = KTrdt(i_k);
    Xrdt = Xrdt + X(:, kTrdt:kTrdt+NTmax-1);
    updateWaitBar(i_k);
end
Xrdt = Xrdt / length(KTrdt);

closeWaitBar();

end

