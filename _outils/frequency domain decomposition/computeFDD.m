function [f, SV, mod_shapes] = computeFDD(t, X, varargin)
%COMPUTEFDD Frequency Domain Decomposition
%   f(freq_index)
%   SV(sv_index, freq_index)
%   mod_shapes(sensor_index, sv_index, freq_index)

p = inputParser;
% paramètres par defaut
NaveragingDef = 1; % moyennage, nombre de découpes du temps initial
MaxLagDef = inf; % troncage de xcorr
XcorrScaleDef = 'unbiased'; % normalisation de l'intercorr, 'biased' ou 'unbiased'
HalfXcorrDef = false; % seulement xcorr avec lag positif

addParameter(p, 'Naveraging', NaveragingDef);
addParameter(p, 'MaxLag', MaxLagDef);
addParameter(p, 'XcorrScale', XcorrScaleDef);
addParameter(p, 'HalfXcorr', HalfXcorrDef);

parse(p, varargin{:});

Naveraging = p.Results.Naveraging;
MaxLag = p.Results.MaxLag;
XcorrScale = p.Results.XcorrScale;
HalfXcorr = p.Results.HalfXcorr;

if Naveraging > 1 && MaxLag < inf
    warning('moyennage ET troncage');
end

%% tests entrees

% taille matrices
if min(size(t)) ~= 1 % t vecteur
    error('erreur taille t');
end
if size(t, 1) > 1 % t vecteur ligne
    t = t.';
end
Nt0 = length(t);

if size(X, 2) ~= Nt0
    X = X.';
end
if size(X, 2) ~= Nt0
    error('erreur taille t ou X');
end
Nx = size(X, 1);

% pas de temps
dt = (t(end)-t(1))/(Nt0-1);

if any(abs(diff(t)/dt-1) > 1e-4)
    error('pas de temps non constant');
end


%% calcul xcorr

% decoupage temporel
Ntav = floor(Nt0/Naveraging); % taille sous-signaux

% max lag
Nmaxlag = round(MaxLag/dt);
Nmaxlag = min(Nmaxlag, Ntav-1);

% calcul xcorr
Mcorr = zeros(Nx^2, 2*Nmaxlag+1);
for kav = 0:Naveraging-1
    Mcorr = Mcorr + xcorr(X(:, kav*Ntav+1:(kav+1)*Ntav).', Nmaxlag, XcorrScale).';
end
Mcorr = Mcorr/Naveraging;

if HalfXcorr % partie temporelle droite seulement
    Mcorr = Mcorr(:, Nmaxlag+1:end);
end

%% calcul fft

FFTcorr = fft(Mcorr.').';

% frequences positives
FFTcorr = FFTcorr(:, 1:floor(end/2)+1);
f = (0:size(FFTcorr, 2)-1) / (size(Mcorr, 2)*dt);

% reecriture
FFTcorr = reshape(FFTcorr, [Nx, Nx, size(FFTcorr, 2)]);


%% calcul SVD

SV = nan(Nx, size(FFTcorr, 3)); % SV(sv_index, freq_index)
mod_shapes = nan(size(FFTcorr)); % mod_shapes(sensor_index, sv_index, freq_index)

for kf = 1:size(FFTcorr, 3)
    [mod_shapes(:, :, kf), S, ~] = svd(FFTcorr(:, :, kf).');
    SV(:, kf) = diag(S);
end


end

