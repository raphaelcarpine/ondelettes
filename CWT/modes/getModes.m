function [t, freqs, shapes] = getModes(Ridges, ddlRef)
%GETMODES Summary of this function goes here
%   Ridges = {ridgeddl1, ...}
%   t = {t1 = [...], ...}
%   amplitudes = {ampl1 = [...], ...}
freqTol = 1e-2;

if nargin < 2
    ddlRef = 1;
end

nddl = length(Ridges);
ridgesRef = Ridges{ddlRef};
nridges = length(ridgesRef.freq);

t = cell(1, nridges);
freqs = cell(1, nridges);
shapes = cell(1, nridges);


for kr = 1:nridges % reference ridge increment
    tr = ridgesRef.time{kr};
    ntr = length(tr);
    freqsr = ridgesRef.freq{kr};
    shapesr = ones(1, ntr);
    shapesRef = ridgesRef.val{kr}; % normalization
    
    freqsr = [nan(ddlRef-1, ntr); freqsr; nan(nddl-ddlRef, ntr)];
    shapesr = [nan(ddlRef-1, ntr); shapesr; nan(nddl-ddlRef, ntr)];
    
    for kddl = [1:(ddlRef-1), (ddlRef+1):nddl] % ddl increment
        ridges2 = Ridges{kddl};
        
        for kr2 = 1:length(ridges2.time) % ridge increment
            if ridges2.time{kr2}(end) < ridgesRef.time{kr}(1) || ridges2.time{kr2}(1) > ridgesRef.time{kr}(end)
                continue
            end
            
            kt2 = 1; % time ind of other ridge
            for kt = 1:ntr % time ind of ref ridge
                while kt2 <= length(ridges2.time{kr2}) && ridges2.time{kr2}(kt2) < ridgesRef.time{kr}(kt)
                    kt2 = kt2+1;
                end
                if kt2 > length(ridges2.time{kr2})
                    continue
                end
                
                if abs(ridges2.freq{kr2}(kt2) - ridgesRef.freq{kr}(kt)) / ridgesRef.freq{kr}(kt) < freqTol
                    freqsr(kddl, kt2) = ridges2.freq{kr2}(kt2);
                    shapesr(kddl, kt2) = ridges2.val{kr2}(kt2);
                end
            end
        end
        
        % normalization
        shapesr(kddl, :) = shapesr(kddl, :) ./ shapesRef;
    end
    
    
    t{kr} = tr;
    freqs{kr} = freqsr;
    shapes{kr} = shapesr;
end

end