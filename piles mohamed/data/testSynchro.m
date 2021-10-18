function [synchro, maxDesync, resync, resyncVect, maxDesyncAfterResync] = testSynchro(T, X, channelNames)
%TESTSYNCHRO Summary of this function goes here
%   Detailed explanation goes here

tolSync = 2; % k, k*\Delta t
threshold = 0.2;
tau = 1.3;

X = X - mean(X, 2);

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


if false % test
    shockIndexesMat = cell2mat(shockIndexes');
    shockIndexesMat(:, 1) = shockIndexesMat(:, 1) - min(shockIndexesMat(:, 1)) + 1;
    for kt = 2:size(shockIndexesMat, 2)
        shockIndexesMat(:, kt) = shockIndexesMat(:, kt) - min(shockIndexesMat(:, kt)) + max(shockIndexesMat(:, kt-1)) + 2;
    end
    figure;
    for kch = 1:size(shockIndexesMat, 1)
        scatter(shockIndexesMat(kch, :), kch*ones(size(shockIndexesMat(kch, :))));
        hold on
    end
    for kt = 1:size(shockIndexesMat, 2)-1
        xline(max(shockIndexesMat(:, kt))+1);
    end
end

%% test synchro

synchro = true;
maxDesync = [];
for kch = 2:length(shockIndexes)
    if length(shockIndexes{kch}) ~= length(shockIndexes{1}) % nombre de chocs différent
        synchro = false;
        maxDesync = inf;
        resync = false;
        resyncVect = [];
        maxDesyncAfterResync = inf;
        return
    end
    if any(abs(shockIndexes{kch} - shockIndexes{1}) > tolSync)
        synchro = false;
    end
    maxDesync(end+1) = max(abs(shockIndexes{kch} - shockIndexes{1}));
end

maxDesync = max(maxDesync);
if isempty(maxDesync)
    maxDesync = nan;
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
    end
end

maxDesyncAfterResync = [];
for kch = 2:length(shockIndexes)
    if length(shockIndexes{kch}) ~= length(shockIndexes{1}) % nombre de chocs différent
        maxDesyncAfterResync = inf;
        break
    end
    maxDesyncAfterResync(end+1) = max(abs(shockIndexes{kch} - shockIndexes{1}));
end

maxDesyncAfterResync = max(maxDesyncAfterResync);
if isempty(maxDesyncAfterResync)
    maxDesyncAfterResync = nan;
end

%% test

if false % test
    shockIndexesMat = cell2mat(shockIndexes');
    shockIndexesMat(:, 1) = shockIndexesMat(:, 1) - min(shockIndexesMat(:, 1)) + 1;
    for kt = 2:size(shockIndexesMat, 2)
        shockIndexesMat(:, kt) = shockIndexesMat(:, kt) - min(shockIndexesMat(:, kt)) + max(shockIndexesMat(:, kt-1)) + 2;
    end
    figure;
    for kch = 1:size(shockIndexesMat, 1)
        scatter(shockIndexesMat(kch, :), kch*ones(size(shockIndexesMat(kch, :))));
        hold on
    end
    for kt = 1:size(shockIndexesMat, 2)-1
        xline(max(shockIndexesMat(:, kt))+1);
    end
end

end

