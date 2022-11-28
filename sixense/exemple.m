%% construction de faux signaux
% On construit ici des faux signaux, d'une structure linéaire à trois
% modes, avec 5 capteurs accélérométriques, excitée par un bruit blanc.

% données simulation
T = 10000; % temps total d'acquisition
fe = 20; % fréquence d'échantillonnage
dt = 1/fe; % pas de temps
t = 0:dt:T; % vecteur temps
N = length(t); % nb de points

% données système
f1 = 2; % fréquences propres
f2 = 5;
f3 = 6;
w1 = 2*pi*f1; w2 = 2*pi*f2; w3 = 2*pi*f3;
z1 = 0.01; % taux d'amortissement
z2 = 0.01;
z3 = 0.02;
phi1 = [1;2;3;2;1]; % déformées modales
phi2 = [1;1;0;-1;-1];
phi3 = [0;0;1;2;1];
A1 = 1; % facteurs de participation modale
A2 = 1;
A3 = 1;
B = randn(1, N); % excitation

% calcul de la réponse
H1 = @(p) A1./(w1^2 + 2*z1*w1*p + p.^2); % fonctions de transfert des modes
H2 = @(p) A2./(w2^2 + 2*z2*w2*p + p.^2);
H3 = @(p) A3./(w3^2 + 2*z3*w3*p + p.^2);
iW = 2i*pi*fe/N*[0:floor(N/2), -ceil(N/2)+1:-1]; % vecteur des pusations
fftB = fft(B); % TF de l'excitation
fftX = phi1*(H1(iW).*fftB) + phi2*(H2(iW).*fftB) + phi3*(H3(iW).*fftB); % TF de la réponse
X = ifft(fftX, [], 2); % position
fftA = iW.^2 .* fftX;
A = ifft(fftA, [], 2); % accélération (signal des capteurs accéléromètres)


%% test

% figure;
% plt = plot(t, A);
% xlabel('Temps [s]');
% ylabel('Accélération [m/s²]');
% WaveletMenu('WaveletPlot', plt);


%% projection sur le mode 1 (2 Hz)
% Construction de la projection sur le mode 1, afin d'éliminer le mode 2.
% On suppose ici que les déformées modales ont déjà été estimées à partir
% des données.

proj = projectionMode(phi1, [phi2, phi3]);
a1 = proj.' * A; % signal accélérométrique ne contenant que le mode 1

% figures fourier
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
MotherWavelet = 'morlet'; % 'morlet', 'cauchy'
ct = 3; % ne pas toucher
Q = 2; % régler à 2, ou plus si mode proche non éliminé
fmin = 1; % f1 - 1
fmax = 3; % f1 + 1
freqs = linspace(fmin, fmax, 300);

% paramètres arete
ridgeContinuity = 'none'; % 'none', 'slope'
slopeTimeConst = 3; % inutile pour ridgeContinuity = 'none'

% calcul
CWT = WvltComp(t, a1, freqs, Q, 'MotherWavelet', MotherWavelet, 'DisplayWaitBar', true); % calcul TOC
ridge = SingleRidgeExtract(t, freqs, CWT, MotherWavelet, Q, ct, ridgeContinuity, slopeTimeConst); % calcul ridge


%% données régressions
% On renomme les variables d'intérêt, et on harmonise la longuer des
% signaux (ceux de l'arête sont amputés des effets de bord)

% données ondelette
t_reg = ridge.time; % "reg" pour regression
freq_reg = ridge.freq; % fréquence instantanée du mode 1
ampl_reg = abs(ridge.val); % amplitude de vibration du mode 1

% données de capteurs de position virtuels
capteur_position1 = X(3, :); % faux capteur de position au ddl 3
capteur_position2 = X(5, :); % faux capteur de position au ddl 5
kti = find(t == t_reg(1)); % indice de temps de début du signal de régression (en dehors des effets de bords)
ktf = find(t == t_reg(end)); % indice de temps de fin du signal
capteur_position1_reg = capteur_position1(kti:ktf);
capteur_position2_reg = capteur_position2(kti:ktf);


%% régressions



















