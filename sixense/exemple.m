%% chargement données
% données de simulation d'un pont sous traffic ambiant, cas sans fissures

load('data_simu.mat'); % Acapt, Ycapt, Temp, fe

A = Acapt; % capteurs virtuels d'accélération, à 1/4, 1/2 et 3/4 de la travée
X = Ycapt; % capteur virtuel de position, à mi-travée

% déformées modales (connues)
phi1 = sin((1:3)*pi/4).';
phi2 = sin((1:3)*pi/2).';

% données complémentaires
N = size(Acapt, 2); % nb de points temporels
dt = 1/fe; % pas de temps
T = N*dt; % temps total d'acquisition
t = dt * (0:N-1); % vecteur temps


%% projection sur le mode 1 (2 Hz)
% Construction de la projection sur le mode 1, afin d'éliminer le mode 2.
% On suppose ici que les déformées modales ont déjà été estimées à partir
% des données.

proj = projectionMode(phi1, phi2); % calcul projection optimale pour isoler le mode 1
a1 = proj.' * A; % signal accélérométrique ne contenant que le mode 1

% figures fourier
fftA = fft(A, [], 2);
ffta1 = fft(a1);
figure;
ax1 = subplot(2,1,1);
plot(fe/N*(0:floor(N/2)), abs(fftA(:, 1:floor(N/2)+1)));
ylabel('Amplitude');
title('Transformées de Fourier des différents capteurs');
ax2 = subplot(2,1,2);
plot(fe/N*(0:floor(N/2)), abs(ffta1(:, 1:floor(N/2)+1)));
xlabel('Fréquence [Hz]');
ylabel('Amplitude');
title('Transformée de Fourier de la projection');
linkaxes([ax1,ax2],'x');


%% calcul de l'arête de la TOC

% paramètres ondelettes
MotherWavelet = 'morlet'; % 'morlet', 'cauchy' (au choix, ne change pas les résultats finaux d'après les simulations)
ct = 3; % ne pas toucher
Q = 2; % régler à 2, ou plus si mode proche non éliminé (à tester, peut changer les résultats finaux, idéal à 2 d'après les simulations)
fmin = 1.08; % f1 - 1
fmax = 3.08; % f1 + 1
freqs = linspace(fmin, fmax, 100); % vecteur des fréquences de la TOC

% paramètres arete
ridgeContinuity = 'none'; % 'none', 'slope' ('none' pour la def par argmax, 'slope' pour la def par optimisation de [Carmona et al., 1997] présentée en annexe G de la thèse)
slopeTimeConst = 3; % temps caractéristique de "lissage" de l'arete, inutile pour ridgeContinuity = 'none'

% calcul arête 5min par 5min
freq_arete = nan(size(t)); % arete (fréquence instantanée)
ampl_arete = nan(size(t)); % squelette (amplitude instantanée)
kti = 1; ktf = 1; % indices début et fin calcul
[~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet); % dispersion temporelle TOC
kt_effets_bord = ceil(fe * ct * DeltaT(fmin)); % largeur effets de bord
[initWaitBar, updateWaitBar, closeWaitBar] = ... % barre de progression
    getWaitBar(N, 'windowTitle', 'Calcul TOC');
initWaitBar(); % barre de progression
while ktf < N
    % calcul indices début et fin heure
    kti = max(ktf - 2*kt_effets_bord, 1); % début à la fin de l'intervalle précédent, avec intervalles de recouvrement d'effets de bord
    ktf = min(ktf + 5*60*fe, N); % fin 5 min plus tard
    
    % calcul TOC et arete
    TOC = WvltComp(t(kti:ktf), a1(kti:ktf), freqs, Q,...
        'MotherWavelet', MotherWavelet, 'DisplayWaitBar', false); % calcul TOC
    arete = SingleRidgeExtract(t(kti:ktf), freqs, TOC, MotherWavelet,...
        Q, ct, ridgeContinuity, slopeTimeConst); % calcul ridge
    
    % enregistrement
    kti_arete = find(t == arete.time(1)); % indice début arête
    ktf_arete = find(t == arete.time(end)); % indice fin arête
    freq_arete(kti_arete:ktf_arete) = arete.freq;
    ampl_arete(kti_arete:ktf_arete) = abs(arete.val); % module, car arete.val est complexe car contient également la phase
    
    updateWaitBar(ktf); % barre de progression
end
closeWaitBar(); % barre de progression


%% décorrélation température

% remplissage données NaN (arête non définie), pour simplification calcul xcorr
freq_arete_nonan = freq_arete;
freq_arete_nonan(isnan(freq_arete)) = mean(freq_arete, 'omitnan'); % on remplace par la moyenne

% temps de déphasage Tau0
Rft = xcorr(freq_arete_nonan - mean(freq_arete_nonan), Temp - mean(Temp), 'biased'); % calcul corrélation
[~, kTau0] = max(abs(Rft)); % argmax
kTau0 = kTau0 - length(freq_arete_nonan); % correction car le déphasage 0 n'est pas à 0 dans le vecteur Rft
if kTau0 < 0
    warning('déphasage température négatif'); % en principe, le déphasage est positif
    kTau0 = 0;
end

% affichage
disp('déphasage température :');
Tau0 = seconds(kTau0*dt);
Tau0.Format = 'hh:mm:ss';
disp(Tau0);
figure;
plot(dt*((1:length(Rft))-length(freq_arete)), Rft);
hold on
plot(dt*kTau0, Rft(kTau0 + length(freq_arete_nonan)), 'r+', 'LineWidth', 1.5);
xlabel('Déphasage [s]');
ylabel('Intercorrélation');

% régression linéaire température
Temp2 = [nan(1, kTau0), Temp(1:end-kTau0)]; % décalage temporel température
I = ~isnan(freq_arete) & ~isnan(Temp2); % sélection indices de tps où l'arête est définie
coeffs = [ones(size(Temp2(I).')), Temp2(I).' - mean(Temp2(I))] \ freq_arete(I).'; % régression affine
betaT = coeffs(2); % coefficient régression
freq_arete_corr = freq_arete - betaT*(Temp2 - mean(Temp2(I))); % fréquence corrigée

% affichage
disp('coefficient corrélation fréquence/température [Hz/°C] :');
disp(betaT);

% % affichage
% figure;
% plot(Temp2(I), freq_arete(I), '+');
% hold on
% plot(get(gca, 'XLim'), coeffs(1) + betaT*(get(gca, 'XLim') - mean(Temp2(I))), 'r--');
% xlabel('Température [°C]');
% ylabel('Fréquence [Hz]');

% affichage
figure;
plot(ampl_arete, freq_arete_corr);
xlabel('Amplitude [m/s²]');
ylabel('Fréquence [Hz]');
figure;
plot(1000*X, freq_arete_corr);
xlabel('Flèche [mm]');
ylabel('Fréquence [Hz]');


%% régressions

% sur tout le signal (24h)
I = ~isnan(freq_arete_corr) & ~isnan(ampl_arete) & ~isnan(X); % sélection indices de tps où l'arête est définie
coeffs = [ones(size(ampl_arete(I).')), ampl_arete(I).', X(:, I).'] \ freq_arete_corr(I).'; % régression affine multi-variables
f0 = coeffs(1);
betaA = coeffs(2);
betaX = coeffs(3:end); % betaX peut être de dimension > 1

% par heure
Nsep = 24; % nb de sous-intervalles de temps
Nint = floor(N/Nsep); % taille intervalle
Coeffs = nan(length(coeffs), Nsep);
for ksep = 1:Nsep
    Isep = 1:N > (ksep-1)*Nint & 1:N <= ksep*Nint &... % sélection indices temps dans l'intervalle
        ~isnan(freq_arete_corr) & ~isnan(ampl_arete) & ~isnan(X); % sélection indices de tps où l'arête est définie
    Coeffs(:, ksep) = [ones(size(ampl_arete(Isep).')), ampl_arete(Isep).', X(:, Isep).']...
        \ freq_arete_corr(Isep).'; % régression affine multi-variables
end
meanCoeffs = mean(Coeffs, 2); % moyennes coefficients
errCoeffs = std(Coeffs, [], 2) / sqrt(Nsep); % bornes erreur coefficient

% affichage
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp('freq = f0 + betaA*A + betaX*X + correction_température + dispersion_stat\,');
disp(' - signal entier :');
disp(['     f0 = ', num2str(f0)]);
disp(['     betaA = ', num2str(betaA)]);
disp(['     betaX = ', num2str(betaX)]);
disp([' - découpage en ', num2str(Nsep), ' sous-signaux :']);
disp(['     f0 = ', num2str(meanCoeffs(1)), ' +- ', num2str(errCoeffs(1))]);
disp(['     betaA = ', num2str(meanCoeffs(2)), ' +- ', num2str(errCoeffs(2))]);
disp(['     betaX = ', num2str(meanCoeffs(3:end)), ' +- ', num2str(errCoeffs(3:end))]);












