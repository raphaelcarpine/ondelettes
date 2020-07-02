function [t, freqs, shapes, amplitudes] = getModesCrossCorr(T, SVrx, SVvectrx, Q, fmin, fmax, NbFreq, Nsv, varargin)
%GETMODES Summary of this function goes here
%   
%   Nsv : nb singular values
%
%   varargin -> ridgeExtract
%
%   t = {{t1 = [...],...},...
%   t{ksvd}{kridge}
%   ...

if size(T, 1) ~= 1
    T = transpose(T);
end

if size(T, 1) ~= 1
    error(['array size problem (', num2str(size(T)), ')']);
end


Ndof = size(SVvectrx{1}, 1);


%%

arrayFreqs = linspace(fmin, fmax, NbFreq);

t = cell(1, Nsv);
freqs = cell(1, Nsv);
shapes = cell(1, Nsv);
amplitudes = cell(1, Nsv);

for ksv = 1:Nsv
    
    %% ridge extract
    
    Ridges = RidgeExtract(T, nan, Q, fmin, fmax, NbFreq, 'Wavelet', SVrx{ksv}, varargin{:});
    
    %%
    Nridges = length(Ridges.time);
    
    t{ksv} = cell(1, Nridges);
    freqs{ksv} = cell(1, Nridges);
    shapes{ksv} = cell(1, Nridges);
    amplitudes{ksv} = cell(1, Nridges);
    
    
    for kr = 1:Nridges
        t{ksv}{kr} = Ridges.time{kr};
        freqs{ksv}{kr} = Ridges.freq{kr};
        amplitudes{ksv}{kr} = Ridges.val{kr};
        
        Nt = length(t{ksv}{kr}); % time indexes
        kt0 = 1;
        while kt0 <= length(T) && T(kt0) < t{ksv}{kr}(1)
            kt0 = kt0+1;
        end
        timeInd = kt0:(kt0 + Nt);
        
        shapes{ksv}{kr} = nan(Ndof, Nt);
        
        for kt = 1:Nt
            kt0 = timeInd(kt); % indice sur l'échelle de temps originale
            
            f = freqs{ksv}{kr}(kt);
            
            if f < fmin
                warning('frequence out of bounds');
                kf1 = 1;
                kf2 = 2;
            elseif f > fmax
                warning('frequence out of bounds');
                kf1 = NbFreq;
                kf2 = NbFreq-1;
            else
                kf2 = 1;
                while arrayFreqs(kf2) < f
                    kf2 = kf2 + 1;
                end
                kf1 = kf2 - 1;
            end
            f1 = arrayFreqs(kf1);
            f2 = arrayFreqs(kf2);
            c1 = (f2-f) / (f2-f1);
            c2 = (f-f1) / (f2-f1);
            
            shape = c1*SVvectrx{ksv}(:, kf1, kt0) + c2*SVvectrx{ksv}(:, kf2, kt0);
            
            ampl = sqrt(transpose(shape) * shape);
            if kt > 1 && norm(-shape-shapes{ksv}{kr}(:, kt-1)) < norm(shape-shapes{ksv}{kr}(:, kt-1)) % shape continuity
                ampl = -ampl;
            end
            
            shape = shape / ampl;
            
            if kt == 1 % shape orientation
                [~, indMax] = max(abs(shape));
                
                indOrientation = 1; %indMax
                
                shape = shape * sign(real(shape(indOrientation)));
            end
            
            shapes{ksv}{kr}(:, kt) = shape;
        end
        
    end
    
end

end














