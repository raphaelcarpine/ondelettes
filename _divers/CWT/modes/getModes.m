function [t, freq, freqs, shapes, amplitudes, errors, ridgesNumber] = getModes(Ridges, ddlRef)
%GETMODES Summary of this function goes here
%   Ridges = {ridgeddl1, ...}
%   t = {t1 = [...], ...}
%   ...

% si true, prend seulement la partie où il y a des ridges pour tous les
% ddl, sauf les ddl qui n'ont pas de ridge du tout pour cette freq
continuite = true;

freqTol = 5e-2;

if nargin < 2
    ddlRef = 1;
end

nddl = length(Ridges);
ridgesRef = Ridges{ddlRef};
nridges = length(ridgesRef.freq);

t = cell(1, nridges);
freq = cell(1, nridges);
freqs = cell(1, nridges);
shapes = cell(1, nridges);
amplitudes = cell(1, nridges);
errors = cell(1, nridges);
ridgesNumber = cell(1, nridges);


for kr = 1:nridges % reference ridge increment
    tr = ridgesRef.time{kr};
    ntr = length(tr);
    freqsr = ridgesRef.freq{kr};
    shapesr = ridgesRef.val{kr}; % normalization
    
    freqsr = [nan(ddlRef-1, ntr); freqsr; nan(nddl-ddlRef, ntr)];
    shapesr = [nan(ddlRef-1, ntr); shapesr; nan(nddl-ddlRef, ntr)];
    errorsr = freqTol * ones(nddl, ntr);
    amplitudesr = nan(1, ntr);
    ridgesNumberr = zeros(nddl, ntr);
    
    ridgesNumberr(ddlRef, :) = kr;
    errorsr(ddlRef, :) = 0;
    
    for kddl = [1:(ddlRef-1), (ddlRef+1):nddl] % ddl increment
        ridges2 = Ridges{kddl};
        
        for kr2 = 1:length(ridges2.time) % ridge increment
            if ridges2.time{kr2}(end) <= ridgesRef.time{kr}(1) || ridges2.time{kr2}(1) >= ridgesRef.time{kr}(end)
                continue
            end
            
            kt2 = 1; % time ind of other ridge
            for kt = 1:ntr % time ind of ref ridge
                while kt2 <= length(ridges2.time{kr2}) && ridges2.time{kr2}(kt2) < ridgesRef.time{kr}(kt)
                    kt2 = kt2+1;
                end
                
                % overlap
                if kt2 > length(ridges2.time{kr2})
                    continue
                end
                if ridges2.time{kr2}(kt2) > ridgesRef.time{kr}(kt)
                    continue
                end
                
                if abs(ridges2.freq{kr2}(kt2) - ridgesRef.freq{kr}(kt)) < errorsr(kddl, kt) * ridgesRef.freq{kr}(kt)
                    errorsr(kddl, kt) = abs(ridges2.freq{kr2}(kt2) - ridgesRef.freq{kr}(kt)) / ridgesRef.freq{kr}(kt);
                    freqsr(kddl, kt) = ridges2.freq{kr2}(kt2);
                    shapesr(kddl, kt) = ridges2.val{kr2}(kt2);
                    ridgesNumberr(kddl, kt) = kr2;
                end
            end
        end
    end
    
    % normalization
    for kt = 1:ntr % time ind of ref ridge
        angleMod = shapesr(ddlRef, kt) / abs(shapesr(ddlRef, kt));
        shapesr(:, kt) = shapesr(:, kt) / angleMod;
        
        shape = shapesr(:, kt);
        for ks = 1:length(shape)
            if isnan(shape(ks))
                shape(ks) = 0;
            end
        end
        amplitudesr(kt) = sqrt( transpose(shape) * shape);
        shapesr(:, kt) = shapesr(:, kt) / amplitudesr(kt);
        
        amplitudesr(kt) = amplitudesr(kt) * angleMod;
    end
    
    % continuite
    if continuite
        % détermination du plus long ridge commun continu
        kt0 = 1;
        T = 0;
        kt0max = kt0;
        Tmax = T;
        for kt = 1:length(tr)
            if all (ridgesNumberr(:, kt) == ridgesNumberr(:, kt0))
                T = T+1;
                if T > Tmax
                    kt0max = kt0;
                    Tmax = T;
                end
            else
                kt0 = kt+1;
                T = 0;
            end
        end
        contRidgeInd = kt0max:(kt0max + Tmax-1);
        
        tr = tr(contRidgeInd);
        freqsr = freqsr(:, contRidgeInd);
        shapesr = shapesr(:, contRidgeInd);
        amplitudesr = amplitudesr(contRidgeInd);
        errorsr = errorsr(:, contRidgeInd);
        ridgesNumberr = ridgesNumberr(:, contRidgeInd);
    end
    
    t{kr} = tr;
    freqs{kr} = freqsr;
    shapes{kr} = shapesr;
    amplitudes{kr} = amplitudesr;
    errors{kr} = errorsr;
    ridgesNumber{kr} = ridgesNumberr;
    
    %mode frequency
    dphi = angle(amplitudesr(2:end) ./ amplitudesr(1:end-1));
    dphi = dphi ./ diff(tr);
    dphi = [dphi(1), (dphi(1:end-1)+dphi(2:end))/2, dphi(end)];
    freq{kr} = dphi/(2*pi);
end

end














