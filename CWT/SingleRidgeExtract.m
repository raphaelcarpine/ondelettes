function ridge = SingleRidgeExtract(t, freqs, CWT, MotherWavelet, Q, ct, ridgeContinuity)
%SINGLERIDGEEXTRACT Summary of this function goes here
%   Detailed explanation goes here


% ridge computation
switch ridgeContinuity
    case 'none'
        [~, Fridge] = max(abs(CWT), [], 1);
        [Fridge, Aridge] = localMax3Points(freqs([Fridge-1; Fridge; Fridge+1]),...
            CWT([Fridge-1; Fridge; Fridge+1] + [1;1;1] * (0:size(CWT, 2)-1)*size(CWT, 1)));
    case 'simple'
        Fridge = nan(1, size(CWT, 2));
        [~, Fridge(1)] = max(abs(CWT(:, 1)));
        localMax = abs(CWT(1:end-1, :)) < abs(CWT(2:end, :));
        localMax = localMax(1:end-1, :) & ~localMax(2:end, :);
        for kt = 2:size(CWT, 2)
            localMaxFreq = find(localMax(:, kt)) + 1;
            [~, closestLocalMax] = min(abs(localMaxFreq - Fridge(kt-1)));
            Fridge(kt) = localMaxFreq(closestLocalMax);
        end
        [Fridge, Aridge] = localMax3Points(freqs([Fridge-1; Fridge; Fridge+1]),...
            CWT([Fridge-1; Fridge; Fridge+1] + [1;1;1] * (0:size(CWT, 2)-1)*size(CWT, 1)));
    otherwise
        error(' ');
end

ridge.time = t;
ridge.val = Aridge;
ridge.freq = Fridge;

% edge effects
[~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
dt = (t(end)-t(1))/(length(t)-1);
Ntedge = ceil(ct*DeltaT(freqs(1))/dt);
nt1 = 1;
for kt = min(Ntedge + 1, length(t)) : -1 : 1
    if ct*DeltaT(ridge.freq(kt)) > (kt-1) * dt
        nt1 = kt;
        break
    end
end
nt2 = length(t);
for kt = min(Ntedge, length(t)-1) : -1 : 0
    if ct*DeltaT(ridge.freq(length(t) - kt)) > kt * dt
        nt2 = length(t)- kt;
        break
    end
end

ridge.time = ridge.time(nt1:nt2);
ridge.val = ridge.val(nt1:nt2);
ridge.freq = ridge.freq(nt1:nt2);

end

