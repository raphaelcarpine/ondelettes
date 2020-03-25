function phase = phaseContinuity(phase)
%PHASECONTINUITY Summary of this function goes here
%   Detailed explanation goes here
for k = 2:length(phase)
    if abs(phase(k) + 2*pi - phase(k-1)) < abs(phase(k) - phase(k-1))
            phase(k:end) = phase(k:end) + 2*pi;
    elseif abs(phase(k) - 2*pi - phase(k-1)) < abs(phase(k) - phase(k-1))
            phase(k:end) = phase(k:end) - 2*pi;
    end
end
end

