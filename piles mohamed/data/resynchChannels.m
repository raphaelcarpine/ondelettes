function [Tshocks, Xshocks, resynched] = resynchChannels(T, X, channelNames)
%RESYNCHCHANNELS0 Summary of this function goes here
%   Detailed explanation goes here

resynch = true;

tolSync = 3; % k, k*\Delta t
threshold = 0.2;
tau = 1.3;

%% sensors channels

sensCh = {};
for kch = 1:length(channelNames)
    if str2double(channelNames{kch}(end-1)) == length(sensCh)
        sensCh{end}(end+1) = kch;
    else
        sensCh{end+1} = kch;
    end
end


%% shock detection

X2 = (X - mean(X, 2)).^2;

sensShocks = cell(size(sensCh));

for ksens = 1:length(sensCh)
    t0 = -inf;
    for kt = 1:size(X, 2)
        if T(kt) < t0 + tau
            continue
        end
        
        if sum(X2(sensCh{ksens}, kt)) >= length(sensCh{ksens}) * threshold^2
            sensShocks{ksens}(end+1) = kt;
            t0 = T(kt);
        end
    end
end


%% shocks matching

sensShocksNb = cellfun(@length, sensShocks);
maxSh = max(sensShocksNb);
for ksens = 1:length(sensShocks)
    sensShocks{ksens} = [nan(1, maxSh), sensShocks{ksens}, nan(1, 2*maxSh - length(sensShocks{ksens}))];
end

for ksens = 2:length(sensShocks)
    while ~isnan(max(abs(sensShocks{ksens} - sensShocks{1}))) % security
        d0 = max(abs(sensShocks{ksens} - sensShocks{1})); % distance
        d1 = max(abs([nan, sensShocks{ksens}(1:end-1)] - sensShocks{1})); % distance decalage 1 vers la droite
        dm1 = max(abs([sensShocks{ksens}(2:end), nan] - sensShocks{1})); % distance decalage 1 vers la gauche
        if (d0 <= d1 || isnan(d1)) && (d0 <= dm1 || isnan(d0))
            break
        elseif d1 < d0
            sensShocks{ksens} = [nan, sensShocks{ksens}(1:end-1)];
        elseif dm1 < d0
            sensShocks{ksens} = [sensShocks{ksens}(2:end), nan];
        end
    end
end

sensShocks = cell2mat(sensShocks.');
while ~isempty(sensShocks) && any(isnan(sensShocks(:, 1)))
    sensShocks = sensShocks(:, 2:end);
end
while ~isempty(sensShocks) && any(isnan(sensShocks(:, end)))
    sensShocks = sensShocks(:, 1:end-1);
end


%% shocks choice

resynched = false;

chosenShocksIndexes = {};
chosenShocksResynch = {};
Nsh = size(sensShocks, 2);
for ksh = 1:Nsh-1
    if max(sensShocks(:, ksh)) - min(sensShocks(:, ksh)) <= tolSync &&... % shock synchronized
            max(sensShocks(:, ksh+1)) - min(sensShocks(:, ksh+1)) <= tolSync &&... % next shock synchronized
            ~any(isnan(X(:, min(sensShocks(:, ksh)):min(sensShocks(:, ksh+1)))), 'all') % no nan values between shocks (time dilatation)
        chosenShocksIndexes{end+1} = [min(sensShocks(:, ksh)), min(sensShocks(:, ksh+1)) - 3];
        chosenShocksResynch{end+1} = zeros(length(sensCh), 1);
        continue
    end
    
    if resynch
        % reference sensor
        [~, refSens] = min(abs(sensShocks(:, ksh) - median(sensShocks(:, ksh))));
        
        % synchronization
        shockResynch = zeros(length(sensCh), 1);
        for ksens = 1:length(sensCh)
            if ksens == refSens
                continue
            end
            
            if abs(sensShocks(ksens, ksh) - sensShocks(refSens, ksh)) > tolSync/2
                shockResynch(ksens) = round(sensShocks(refSens, ksh) -sensShocks(ksens, ksh));
            end
        end
        
        % synch test
        if max(sensShocks(:, ksh) + shockResynch) - min(sensShocks(:, ksh) + shockResynch) <= tolSync &&... % shock synchronized
                max(sensShocks(:, ksh+1) + shockResynch) - min(sensShocks(:, ksh+1) + shockResynch) <= tolSync &&... % next shock synchronized
                ~any(isnan(X(:, min(sensShocks(:, ksh)):min(sensShocks(:, ksh+1)))), 'all') % no nan values between shocks (time dilatation)
            chosenShocksIndexes{end+1} = [min(sensShocks(:, ksh) + shockResynch), min(sensShocks(:, ksh+1) + shockResynch) - 3];
            chosenShocksResynch{end+1} = shockResynch;
            resynched = true;
            continue
        end
    end
end

%% shocks construction

Tshocks = cell(size(chosenShocksIndexes));
Xshocks = cell(size(chosenShocksIndexes));

for ksh = 1:length(Tshocks)
    I = chosenShocksIndexes{ksh};
    Tshocks{ksh} = T(I(1):I(2));
    Xshocks{ksh} = nan(size(X, 1), I(2) - I(1) + 1);
    for ksens = 1:length(sensCh)
        Xshocks{ksh}(sensCh{ksens}, :) = X(sensCh{ksens}, (I(1):I(2)) - chosenShocksResynch{ksh}(ksens));
    end
end


end

