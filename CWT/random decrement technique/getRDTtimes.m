function KTrdt = getRDTtimes(t, X, Tmax, referenceChannel, detectionMethod, detectionThreshold)
%GETRDT Random Decrement Technique
%   detectionMethod = 'threshold', 'zero level positive slope', 'zero level negative slope'

if nargin < 6
    detectionThreshold = 0;
end


%%
dt = mean(diff(t));
NTmax = floor(Tmax / dt) + 1;

switch detectionMethod
    case 'threshold'
        detectNegSlope = true;
        detectPosSlope = true;
    case 'zero level negative slope'
        detectNegSlope = true;
        detectPosSlope = false;
        detectionThreshold = 0;
    case 'zero level positive slope'
        detectNegSlope = false;
        detectPosSlope = true;
        detectionThreshold = 0;
    otherwise
        error('unknown detection method');
end

%% input array size

% array size
if size(t, 1) == 1 && size(t, 2) == size(X, 2)
    dt = mean(diff(t));
elseif size(t, 1) == size(X, 1) && size(t, 2) == size(X, 2)
    dt = mean(diff(t(1, :)));
    for k_t = 2:size(t, 1)
        if max(abs(t(k_t, :) - t(1, :))) > dt * 1e-3
            error(' ');
        end
    end
    t = t(1, :);
else
    error(['array size problem (', num2str(size(t)), ' & ', num2str(size(X)), ')']);
end

% time step
if any(abs(diff(t)/dt - 1) > 1e-3)
    error('non-constant time step');
end


%% times RDT

KTrdt = [];
for kt = 1:length(t)-max(NTmax, 2)+1
    if detectNegSlope && X(referenceChannel, kt) > detectionThreshold...
            && X(referenceChannel, kt+1) <= detectionThreshold
        KTrdt(end+1) = kt;
    end
    if detectPosSlope && X(referenceChannel, kt) <= detectionThreshold...
            && X(referenceChannel, kt+1) > detectionThreshold...
            && kt+1 <= length(t)-NTmax+1
        KTrdt(end+1) = kt+1;
    end
end

end

