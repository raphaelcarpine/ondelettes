%% données pont

L = 17.5; % longueur du pont
E = 30e9; % module d'young du béton
rho = 2500; % masse volumique du béton
l_pont = 4.85; % largeur du pont
h_pont = 1.2;
J = l_pont * h_pont^3 / 12; % moment quadratique du pont
mu = l_pont * h_pont * rho; % masse linéique du pont
amort = 4e5; % coefficient d'amortissement


%% integration spatiale

N = 4+3:350; % nb de points ds l'espace

nbFreqs = 3;


%% calcul freq propre

% frequences therorique, problème continu
freqsTh = pi/(2*L^2) * sqrt(E*J/mu) * (1:nbFreqs).^2;
freqsTh = transpose(freqsTh);

% frequences propres, problème discrétisé
Freqs = nan(nbFreqs, length(N));

for kn = 1:length(N)
    n = N(kn);
    
    % pas spatial
    dx = L / (n-1);
    
    % masse
    M = mu * diag([0, 0, ones(1, n-4), 0, 0]);
    
    % rigidité
    K_eq = zeros(n); % equa diff
    for i = 3:(n-2)
        K_eq(i, i-2:i+2) = [1 -4 6 -4 1];
    end
    K_eq = E*J/dx^4 * K_eq;
    
    K_CL1 = zeros(n); % CL \nu(x=0) = \nu(x=L) = 0
    K_CL1(1, 1) = 1;
    K_CL1(n, n) = 1;
    K_CL1 = 1/dx^4 * K_CL1;
    
    K_CL2 = zeros(n); % CL \partial^3\nu/\partial x^3(x=0) = \partial^3\nu/\partial x^3(x=L) = 0
    K_CL2(2, 1:3) = [-1 2 -1];
    K_CL2(n-1, n-2:n) = [-1 2 -1];
    K_CL2 = 1/dx^4 * K_CL2;
    
    K = K_eq + K_CL1 + K_CL2;
    
    % restriction aux ddl 3 à n-1 (sinon M non inversible)
    % (y1 = yn = 0, et y2 = y3/2 et yn-1 = yn-2/2)
    Mtransition = diag([0, 0, ones(1, n-4), 0, 0]);
    Mtransition(2, 3) = 1/2;
    Mtransition(n-1, n-2) = 1/2;
    
    % matrices réduites
    Mr = M * Mtransition;
    Kr = K * Mtransition;
    Mr = Mr(3:n-2, 3:n-2);
    Kr = Kr(3:n-2, 3:n-2);
    
    %% calcul freq propre, freq excitation
    
    % frequence propre, problème discrétisé
    freqs = eig(Mr\Kr);
    freqs = sqrt(freqs)/(2*pi);
    freqs = sort(freqs);
    freqs = freqs(1:nbFreqs);
    
    Freqs(:, kn) = freqs;
    
end

%% affichage

figure;
plot(N, N .* (Freqs ./ freqsTh - 1));
ylim0 = get(gca, 'YLim');
set(gca, 'YLim', [0, ylim0(2)]);
xlabel('N');
ylabel('N * erreur freq. propre');
set(gca, 'XGrid', 'on', 'YGrid', 'on');

figure;
plot(N, Freqs ./ freqsTh - 1);
xlabel('N');
ylabel('erreur freq. propre');
set(gca, 'YScale', 'log');
set(gca, 'XGrid', 'on', 'YGrid', 'on');

figure;
plot(N, Freqs);
hold on
plot(N, freqsTh * ones(1, length(N)), 'r--');
xlabel('N');
ylabel('freq. propre');

