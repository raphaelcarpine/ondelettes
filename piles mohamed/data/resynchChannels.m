function [T, X] = resynchChannels(T, X, resyncVect)
%RESYNCHCHANNELS Summary of this function goes here
%   Detailed explanation goes here

if size(X, 1) > size(X, 2)
    X = X.';
end

for kch = 1:size(X, 1)
    
end


end

