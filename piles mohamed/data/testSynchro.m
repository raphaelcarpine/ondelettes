function [synchro, resync, resyncVect] = testSynchro(T, X, channelNames)
%TESTSYNCHRO Summary of this function goes here
%   Detailed explanation goes here

tolSync = 1; % k, k*\Delta t
threshold = 0.1;
tau = 1.3;

if size(X, 1) > size(X, 2)
    X = X.';
end

%% shock indexes

shockIndexes = cell(1, size(X, 1));
for kch = 1:size(X, 1)
    t0 = -inf;
    for kt = 1:size(X, 2)
        if T(kt) < t0 + tau
            continue
        end
        
        if abs(X(kch, kt)) >= threshold
            shockIndexes{kch}(end+1) = kt;
            t0 = T(kt);
        end
    end
end

%% test synchro

synchro = true;
for kch = 2:length(shockIndexes)
    if length(shockIndexes{kch}) ~= length(shockIndexes{1}) % nombre de chocs différent
        synchro = false;
        resync = false;
        resyncVect = [];
        return
    end
    if any(abs(shockIndexes{kch} - shockIndexes{1}) > tolSync)
        synchro = false;
        break
    end
end

if synchro
    resync = false;
    resyncVect = [];
    return
end

%% resynch

resync = true;
resyncVect = [];
resyncVect(1) = 0;
for kch = 2:length(shockIndexes)
    resyncVect(kch) = round(mean(shockIndexes{1} - shockIndexes{kch}));
    shockIndexes{kch} = shockIndexes{kch} + resyncVect(kch);
    if any(abs(shockIndexes{kch} - shockIndexes{1}) > tolSync)
        resync = false;
        return
    end
end


end

