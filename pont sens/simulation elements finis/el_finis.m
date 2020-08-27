clear all

save_results = false;

%% données pont

L = 17.5; % longueur du pont
E = 30e9; % module d'young du béton
rho = 2500; % masse volumique du béton
l_pont = 4.85; % largeur du pont
h_pont = 1;
J = l_pont * h_pont^3 / 12; % moment quadratique du pont
mu = l_pont * h_pont * rho; % masse linéique du pont
amort = 4e5; % coefficient d'amortissement

% frequence propre 1
disp(['freq propre 1 : ', num2str(pi/(2*L^2) * sqrt(E*J/mu))]);


%% integration spatiale

N = 49; % nb de points ds l'espace
dx = L / (N-1);

%% données capteurs ponts

pos_capteurs = [L/6, L/3, L/2, 2*L/3, 5*L/6];
pos_capteurs = linspace(0, L, 100); pos_capteurs = pos_capteurs(2:end-1);


%% données train

mu_t = 0.5 * mu; % masse linéique du train
L_wagons = 14; % écart entre les wagons
N_bogies = 1; % nombre de wagons
c = 78.5; % vitesse du train
g = 9.81;

%test
mu = 1e-3 * mu;


%% données temps

ti = -2;
tf = 10;
% t0 = 0;
dt = 0.01;

t = ti:dt:tf;


%% matrices

% masse
M = mu * diag([0, 0, ones(1, N-4), 0, 0]);

% rigidité
K_eq = zeros(N); % equa diff
for i = 3:(N-2)
    K_eq(i, i-2:i+2) = [1 -4 6 -4 1];
end
K_eq = E*J/dx^4 * K_eq;

K_CL1 = zeros(N); % CL \nu(x=0) = \nu(x=L) = 0
K_CL1(1, 1) = 1;
K_CL1(N, N) = 1;
K_CL1 = 1/dx^4 * K_CL1;

K_CL2 = zeros(N); % CL \partial^3\nu/\partial x^3(x=0) = \partial^3\nu/\partial x^3(x=L) = 0
K_CL2(2, 1:3) = [-1 2 -1];
K_CL2(N-1, N-2:N) = [-1 2 -1];
K_CL2 = 1/dx^4 * K_CL2;

K = K_eq + K_CL1 + K_CL2;

% amortissement
C = zeros(N);
for i = 3:N-2
    C(i, i-1:i+1) = [-1 2 -1];
end
C = amort/dx^2 * C;

% restriction aux ddl 3 à N-1 (sinon M non inversible)
% (y1 = yN = 0, et y2 = y3/2 et yN-1 = yN-2/2)
Mtransition = diag([0, 0, ones(1, N-4), 0, 0]);
Mtransition(2, 3) = 1/2;
Mtransition(N-1, N-2) = 1/2;

% matrices réduites
Mr = M * Mtransition;
Cr = C * Mtransition;
Kr = K * Mtransition;
Mr = Mr(3:N-2, 3:N-2);
Cr = Cr(3:N-2, 3:N-2);
Kr = Kr(3:N-2, 3:N-2);


%% influence du train

% % position des essieux
% essieux = L_bogies * (0:N_wagons-1);
% essieux = [essieux, essieux + l_bogies];
% essieux = essieux - max(essieux);
% essieux = sort(essieux);

% position des essieux
essieux = L_wagons * (0:N_bogies-1);
essieux = essieux - max(essieux);
essieux = sort(essieux);

%test
essieux = L/2;
c = 0;

if any(diff(essieux) <= 2*dx)
    warning('écart d''essieux < 2dx');
end

% masse ajoutée
M_ajout = @(t) mu_t*L_wagons/dx * diag( influence_train(essieux + c*t, dx, N));

% excitation pesanteur
F = @(t) -g * mu_t*L_wagons/dx  * transpose( influence_train(essieux + c*t, dx, N));


%% integration numerique

D = @(t, YdY) [ YdY(N-3:end);
    (Mr + M_ajout(t)) \ (-Kr * YdY(1:N-4) - Cr * YdY(N-3:end) + F(t))];

Y0 = zeros(N-4, 1);
V0 = zeros(N-4, 1);
YV0 = [Y0; V0];


options = odeset('OutputFcn', @waitbarOutput);
[t, YV] = ode45(D, t, YV0, options);

t = t';
YV = YV';
Y = YV(1:N-4, :);

Ytot = getYtot(Y);


%% affichage

% animation
moving_coeff = 1.;
movingPlot(Ytot, t, L, essieux, c, pos_capteurs, moving_coeff);

% wavelet capteurs
Ycapt = getYcapt(Ytot, pos_capteurs, dx);
Ycapt = Ycapt + 1e-6 * randn(size(Ycapt));


fig = figure;
ax = axes(fig);
plt = plot(ax, t, Ycapt);

fmin = 3;
fmax = 16;
Q = 10;
MaxRidges = 1;
RealShapePlot = deformeePont(L, pos_capteurs);
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', true,...
    'MaxRidges', MaxRidges, 'RealShapePlot', RealShapePlot);


%% enregistrement

if save_results
    listing = dir('pont sens/simulation elements finis/resultats');
    listingNames = {listing.name};
    nbSimul = 0;
    for kname = length(listingNames)
        name1 = strsplit(listingNames{kname}, '_');
        name1 = name1{1};
        if length(name1) > 5 && strcmp(name1(1:5), 'simul')
            nbSimul = max(nbSimul, str2double(name1(6:end)));
        end
    end
    save(sprintf('pont sens/simulation elements finis/resultats/simul%d_N%d', [nbSimul+1, N]),...
        'L', 'E', 'rho', 'l_pont', 'h_pont', 'J', 'mu', 'amort', 'N', 'dx', 'pos_capteurs', 'mu_t',...
        'L_wagons', 'N_bogies', 'c', 'g', 'ti', 'tf', 'dt', 't', 'Ytot', 'Y0', 'V0',...
        'essieux');
end


%% influence des essieux
function influence = influence_train(essieux_t, dx, N)
essieux_ind = essieux_t / dx + 1;
essieux_ind = essieux_ind(essieux_ind >= 1 & essieux_ind < N);
influence_gauche = zeros(1, N);
influence_gauche(floor(essieux_ind)) = 1 - essieux_ind + floor(essieux_ind);
influence_droite = zeros(1, N);
influence_droite(floor(essieux_ind)+1) = essieux_ind - floor(essieux_ind);
influence = influence_gauche + influence_droite;
influence = influence(3:N-2); % supression des ddl des CL
end



%% interpolation etc

function Ytot = getYtot(Y) % rajoute les DDL 1, 2, N-1 et N
Ytot = [zeros(1, size(Y, 2)); Y(1, :)/2; Y; Y(end, :)/2; zeros(1, size(Y, 2))];
end


%% barre de chargement

function status = waitbarOutput(t, ~, flag)
persistent ti tf f t_last t_1s T_1s

status = 0;
if nargin < 3 || isempty(flag)
    t_0s = 24*3600*now;
    if t_0s - t_last > 0.1
        t_last = t_0s;
        x = (t(end)-ti)/(tf-ti);
        try
            waitbar(x, f, sprintf('t = %.2f', t(end)));
        catch
            status = 1;
        end
    end
    if t_0s - t_1s > 2
        t_remaining = (t_0s - t_1s) * (tf-t(end))/(t(end)-T_1s);
        set(f, 'Name', sprintf('Computing ODE (%dmin remaining)', round(t_remaining/60)));
        t_1s = t_0s;
        T_1s = t(end);
    end
elseif strcmp(flag, 'init')
    f = waitbar(0, sprintf('t = %.2f', t(1)), 'Name', 'Computing ODE');
    ti = t(1);
    tf = t(end);
    t_last = 24*3600*now;
    t_1s = t_last;
    T_1s = t(1);
elseif strcmp(flag, 'done')
    try
        close(f);
    catch
    end
end
end


%%

% K = 
% E*J*
% [1  0  0  0  0  0  0... ; -> CL \nu(x=0) = 0
%  1 -3  3 -1  0  0  0... ; -> CL \partial^3\nu/\partail x^3(x=0) = 0
%  1 -4  6 -4  1  0  0... ; -> eq diff
%  0  1 -4  6 -4  1  0... ; -> eq diff
% ...
%     ... 1 -4  6 -4  1  0  0 ; -> eq diff
%     ... 0  1 -4  6 -4  1  0 ; -> eq diff
%     ... 0  0  1 -4  6 -4  1 ; -> eq diff
%     ... 0  0  0  1 -3  3 -1 ;  -> CL \partial^3\nu/\partail x^3(x=L) = 0
%     ... 0  0  0  0  0  0  1 ;]  -> CL \nu(x=L) = 0
% 
% 
% 
% M = mu*dx*
% [0  0  0  0  0  0  0... ; -> CL \nu(x=0) = 0
%  0  0  0  0  0  0  0... ; -> CL \partial^3\nu/\partail x^3(x=0) = 0
%  0  0  1  0  0  0  0... ; -> eq diff
%  0  0  0  1  0  0  0... ; -> eq diff
% ...
%     ... 0  0  1  0  0  0  0 ; -> eq diff
%     ... 0  0  0  1  0  0  0 ; -> eq diff
%     ... 0  0  0  0  1  0  0 ; -> eq diff
%     ... 0  0  0  0  0  0  0 ;  -> CL \partial^3\nu/\partail x^3(x=L) = 0
%     ... 0  0  0  0  0  0  0 ;]  -> CL \nu(x=L) = 0
% 
% 
% 
% Ensuite, il faut laisser à 0 tous les termes en position 1, 2, N-2 et N-1 dans la matrice C, et dans l'excitation.
% Le plus simple je pense, c'est d'utiliser une matrice, Mddl, qui va enlever ces ddl.
% 
% Mddl = diag([0, 0, ones(1, N-4), 0, 0]);
% M = mu*Mddl*eye(N)
% q = Mddl*q; (q est un vecteur colone)
% mu_q = Mddl*mu_q;  (mu_q est une matrice)
% 
% 
% dis moi si ça marche
% 
% hésite pas si t'as des questions, bon courage
% 
