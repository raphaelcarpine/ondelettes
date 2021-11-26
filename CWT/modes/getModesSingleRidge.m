function [t, freq, shapes, amplitudes] = getModesSingleRidge(X, Y, Q, fmin, fmax, NbFreq, varargin)
%GETMODES Summary of this function goes here
%   X = t
%   Y(ddl, kt)
%
%   t = {t1 = [...], ...}
%   ...

%%

multiSignalModeAbsValueDefault = false;
multiSignalModeAbsValueModeAmplDefault = false;

p = inputParser;
p.KeepUnmatched = true;
p.addOptional('MultiSignalModeAbsValue', multiSignalModeAbsValueDefault); % ridge detection sum|CWT|² instead of |sum CWT²|
p.addOptional('MultiSignalModeAbsValueModeAmpl', multiSignalModeAbsValueModeAmplDefault); % mode amplitude as sqrt(sum|CWT(fr)|²)
p.addOptional('MotherWavelet', '');
p.addOptional('DerivationOrder', nan);

parse(p, varargin{:});

multiSignalModeAbsValue = p.Results.MultiSignalModeAbsValue;
multiSignalModeAbsValueModeAmpl = p.Results.MultiSignalModeAbsValueModeAmpl;
MotherWavelet = p.Results.MotherWavelet;
DerivationOrder = p.Results.DerivationOrder;

varargin2 = namedargs2cell(p.Unmatched);

%%

if size(X, 1) ~= 1
    X = transpose(X);
end
if size(X, 2) ~= size(Y, 2)
    Y = transpose(Y);
end

if size(X, 1) ~= 1 || size(X, 2) ~= size(Y, 2)
    error(['array size problem (', num2str(size(X)), ' & ', num2str(size(Y)), ')']);
end


Nddl = size(Y, 1);

%%
wvlt = cell(1, Nddl);
for kddl = 1:Nddl
    wvlt{kddl} = WvltComp(X, Y(kddl, :), linspace(fmin, fmax, NbFreq), Q,...
        'MotherWavelet', MotherWavelet, 'DerivationOrder', DerivationOrder);
end

wvlt2 = zeros(size(wvlt{1}));
for kddl = 1:Nddl
    if ~multiSignalModeAbsValue
        wvlt2 = wvlt2 + wvlt{kddl}.^2;
    else
        wvlt2 = wvlt2 + abs(wvlt{kddl}).^2;
    end
end

Ridges = RidgeExtract(X, nan, Q, fmin, fmax, NbFreq, 'Wavelet', wvlt2,...
    'MotherWavelet', MotherWavelet, 'SquaredCWT', true, varargin2{:});

%%
Nridges = length(Ridges.time);

t = cell(1, Nridges);
freq = cell(1, Nridges);
shapes = cell(1, Nridges);
amplitudes = cell(1, Nridges);

arrayFreqs = linspace(fmin, fmax, NbFreq);


for kr = 1:Nridges    
    t{kr} = Ridges.time{kr};
    freq{kr} = Ridges.freq{kr};
    
    Nt = length(t{kr});
    kt0 = 1;
    while kt0 <= length(X) && X(kt0) < t{kr}(1)
        kt0 = kt0+1;
    end
    timeInd = kt0:(kt0 + Nt);
    
    shapes{kr} = nan(Nddl, Nt);
    amplitudes{kr} = nan(1, Nt);
    
    for kt = 1:Nt
        kt0 = timeInd(kt); % indice sur l'échelle de temps originale
        
        f = freq{kr}(kt);
        
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
        
        shape = nan(Nddl, 1);
        for kddl = 1:Nddl
            shape(kddl) = c1*wvlt{kddl}(kf1, kt0) + c2*wvlt{kddl}(kf2, kt0);
        end
        
        ampl = sqrt(transpose(shape) * shape);
        
        if kt > 1 && abs(-ampl-amplitudes{kr}(kt-1)) < abs(ampl-amplitudes{kr}(kt-1)) % phase continuity
            ampl = -ampl;
        end
        
        shape = shape / ampl;
        
        if kt == 1 % shape orientation
            [~, indMax] = max(abs(shape));
            
            indOrientation = 1; %indMax
            
            ampl = ampl * sign(real(shape(indOrientation)));
            shape = shape * sign(real(shape(indOrientation)));
        end
        
        shapes{kr}(:, kt) = shape;
        amplitudes{kr}(kt) = ampl;
    end
    
    if multiSignalModeAbsValueModeAmpl
        amplitudes{kr} = sqrt(Ridges.val{kr});
    end
    
end


end














