function [T, X] = resynchChannels2(T, X, resyncVect)
%RESYNCHCHANNELS Summary of this function goes here
%   Detailed explanation goes here

if size(X, 1) > size(X, 2)
    X = X.';
end

for kch = 1:size(X, 1)
    if resyncVect(kch) > 0
        X(kch, 1+resyncVect(kch):end) = X(kch, 1:end-resyncVect(kch));
    elseif resyncVect(kch) < 0
        X(kch, 1:end+resyncVect(kch)) = X(kch, 1-resyncVect(kch):end);
    end
end

X = X(:, 1+max(resyncVect):end+min(resyncVect));
T = T(1+max(resyncVect):end+min(resyncVect));

end

