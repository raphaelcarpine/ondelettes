function [t, freq, freqs, shapes, amplitudes, errors, ridgesNumber] = getModes(Ridges, ddlRef)
%GETMODES Summary of this function goes here
%   Ridges = {ridgeddl1, ...}
%   t = {t1 = [...], ...}
%   ...
freqTol = 1e-2;

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
    ridgesNumberr = nan(nddl, ntr);
    
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
                if kt2 > length(ridges2.time{kr2}) % ridge2 plus court que le ridge de ref à droite
                    continue
                end
                
                if abs(ridges2.freq{kr2}(kt2) - ridgesRef.freq{kr}(kt)) / ridgesRef.freq{kr}(kt) < errorsr(kddl, kt)
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
        amplitudesr(kt) = sqrt(shape.' * shape);
        shapesr(:, kt) = shapesr(:, kt) / amplitudesr(kt);
        
        amplitudesr(kt) = amplitudesr(kt) * angleMod;
    end
    
    %mode frequency
    dphi = angle(amplitudesr(2:end) ./ amplitudesr(1:end-1));
    dphi = dphi ./ diff(tr);
    dphi = [dphi(1), (dphi(1:end-1)+dphi(2:end))/2, dphi(end)];
    freq{kr} = dphi/(2*pi);
    
    t{kr} = tr;
    freqs{kr} = freqsr;
    shapes{kr} = shapesr;
    amplitudes{kr} = amplitudesr;
    errors{kr} = errorsr;
    ridgesNumber{kr} = ridgesNumberr;
end

end














