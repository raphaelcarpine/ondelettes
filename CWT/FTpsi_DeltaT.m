function [FTpsi, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet)
%MU_PSI transformee de fourier de psi et \Delta f
%   FTpsi est normalisé tq max|FTpsi|=1 et argmax|FTpsi|=1
%   DeltaT a pour paramètre f=f0/a=1/a (f0=argmax|FTpsi|)

if strcmp(MotherWavelet, 'cauchy')
    n_cau = (2*Q.^2 - 1/2); % parametre de l'ondelette de Cauchy
    DeltaT_psi = (n_cau/(2*pi)) / sqrt(2*n_cau-1);
    
    FTpsi = @(nu) exp(n_cau*(1 - max(nu, 0)) + n_cau*log(max(nu, 0))) .* (nu > 0); % Calcul Ondelettes cauchy en a * f pour f>0
    DeltaT = @(f) DeltaT_psi ./ f;
    
elseif strcmp(MotherWavelet, 'morlet')
    delta_mor = sqrt(2) * Q / (2*pi); % parametre de l'ondelette de Morlet
    DeltaT_psi = delta_mor / sqrt(2);
    
    FTpsi = @(nu) exp( - (2*pi)^2*(nu-1).^2 * delta_mor^2/2); % Calcul Ondelettes Morlet en a * f
    DeltaT = @(f) DeltaT_psi ./ f;
    
elseif strcmp(MotherWavelet, 'harmonic')
    deltaf_lil = sqrt(3)/Q; % parametre de l'ondelette de harmonique
    DeltaT_psi = 0;
    
    FTpsi = @(nu) (nu >= 1 - deltaf_lil/2) & (nu <= 1 + deltaf_lil/2); % Calcul Ondelettes harmonique en a * f pour f>0
    DeltaT = @(f) DeltaT_psi ./ f;
    
elseif strcmp(MotherWavelet, 'littlewood-paley')
    deltaf_lil = sqrt(3)/Q; % parametre de l'ondelette de Littlewood-Paley
    DeltaT_psi = 0;
    
    FTpsi = @(nu) ((nu >= 1 - deltaf_lil/2) & (nu <= 1 + deltaf_lil/2)) |...
        ((nu >= -1 - deltaf_lil/2) & (nu <= -1 + deltaf_lil/2)); % Calcul Ondelettes Littlewood-Paley en a * f
    DeltaT = @(f) DeltaT_psi ./ f;
    
end
end

