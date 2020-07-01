function [t, freqs, shapes, amplitudes] = getModesCrossCorr(T, X, Q, fmin, fmax, NbFreq, Nsv, varargin)
%GETMODES Summary of this function goes here
%   X = t
%   Y(ddl, kt)
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
if size(T, 2) ~= size(X, 2)
    X = transpose(X);
end

if size(T, 1) ~= 1 || size(T, 2) ~= size(X, 2)
    error(['array size problem (', num2str(size(T)), ' & ', num2str(size(X)), ')']);
end


Ndof = size(X, 1);

%% high pass filtering

% TODO

%% cross correlation & wvlt & SVD

Rx = crossCorrelation(X);

[SVrx, SVvectrx] = svdCWT(T, Rx, fmin, fmax, NbFreq, Q, Nsv);

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
    
    t{sv} = cell(1, Nridges);
    freqs{sv} = cell(1, Nridges);
    shapes{sv} = cell(1, Nridges);
    amplitudes{sv} = cell(1, Nridges);
    
    
    for kr = 1:Nridges
        t{sv}{kr} = Ridges.time{kr};
        freqs{sv}{kr} = Ridges.freq{kr};
        amplitudes{sv}{kr} = Ridges.val{kr};
        
        Nt = length(t{sv}{kr}); % time indexes
        kt0 = 1;
        while kt0 <= length(X) && X(kt0) < t{sv}{kr}(1)
            kt0 = kt0+1;
        end
        timeInd = kt0:(kt0 + Nt);
        
        shapes{sv}{kr} = nan(Ndof, Nt);
        
        for kt = 1:Nt
            kt0 = timeInd(kt); % indice sur l'échelle de temps originale
            
            f = freqs{kr}(kt);
            
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
            if kt > 1 && norm(-shape-shapes{ksv}{kr}(:, kt-1)) < norm(-shape-shapes{ksv}{kr}(:, kt-1)) % shape continuity
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














