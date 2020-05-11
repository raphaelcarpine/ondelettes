
%%%%LOAD Workspace%%%%
load('modes_avant.mat')


figure
plt_p07_avant_mean = plot(temps2_avant, a_p07_avant_mean);
grid on
xlabel('Temps [s]');
ylabel('Accélération Verticale')
title('Mesures Accéléromètre Capteur A7 - Pont de Sens - Avant - 24/06/2003')

% ondelette - TO
WaveletMenu('WaveletPlot', plt_p07_avant_mean, 'fmin', 1, 'fmax', 10, 'Q', 10, 'MaxRidges', 5,'MaxParallelRidges', 5, 'WvltAxesTitle', 'T.O. - Pont de Sens - Avant - 24/06/2003');

%%%%%%%%%%%Modes

%% CWT
NbMaxRidges = 2;
NbMaxParallelRidges = inf;
fmin = 6;
fmax = 20;
nbFreqs = 300;
Q = 8.3;
t = temps2_avant;
X = a_p07_avant_mean;
X_tr = transpose(a_p07_avant_mean);
%%%Peut-etre regarder s'il faudrait prendre le vecteur d'acceleration a_p07
%%%transpose -- Attention car cela a crashé mon matlab!

[time, freq, shapes, amplitudes] = getModesSingleRidge(t, X, Q, fmin, fmax, nbFreqs,... % attention t et time sont différents ! time est le vecteur de temps du ridge uniquement
    'NbMaxRidges', NbMaxRidges, 'NbMaxParallelRidges', NbMaxParallelRidges);
