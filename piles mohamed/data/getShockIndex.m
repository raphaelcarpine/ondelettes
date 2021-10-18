function Kt0 = getShockIndex(T, X, channelNames)
%GETSHOCKINDEX Summary of this function goes here
%   Detailed explanation goes here

refChannels = {'2x', '2y'};
threshold = 1; % m/s²
tau = 1.5; % s

X = X - mean(X, 2);

%%

threshold = length(refChannels) * threshold^2;

I = false(1, length(channelNames));
for kch = 1:length(channelNames)
    if any(strcmp(refChannels, channelNames{kch}(end-1:end)))
        I(kch) = true;
    end
end
X = X(I, :);
X = sum(X.^2, 1);

%%

Kt0 = [];
t0 = -inf;
for kt = 1:length(X)
    if T(kt) < t0 + tau
        continue
    end
    
    if X(kt) >= threshold
        Kt0(end+1) = kt;
        t0 = T(kt);
    end
end
end

