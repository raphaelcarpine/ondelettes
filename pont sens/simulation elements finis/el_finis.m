%% données pont

L = 17.5; % longueur du pont
L = 10;
E = 1; % module d'young du béton
J = 1; % moment quadratique du pont
mu = 1; % masse linéique du pont
amort = 1;


%% données train

mu_t = 1; % masse linéique du train
L_bogies = 18.5; % écartement entre les bogies
l_bogies = 1; % écartement des deux essieux au sein d'un bogie
N_wagons = 16; % nombre de wagons
c_train = 100; % vitesse du train


%% integration spatiale

N = 11; % nb de points ds l'espace
dx = L / (N-1);
%x = linspace(0, L, N);


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

% position des essieux
essieux = L_bogies * (0:N_wagons-1);
essieux = [essieux, essieux + l_bogies];
essieux = essieux - max(essieux);
essieux = sort(essieux);

% masse ajoutée
M_ajout = @(t) mu_t * diag( influence_train(essieux - c*t, dx, N));

% excitation pesanteur
F = @(t) -g*mu_t * transpose( influence_train(essieux - c*t, dx, N));


%% integration numerique

D = @(t, YdY) [ YdY(N-3:end);
    (Mr + M_ajout(t)) \ (-Kr * YdY(1:N-4) - Cr * YdY(N-3:end) + F(t))];




%% influence des essieux
function influence = influence_train(essieux_t, dx, N)
essieux_ind = essieux_t / dx + 1;
essieux_ind = essieux_ind(essieux_ind >= 1 & essieux_ind < N);
influence_gauche = zeros(1, N);
influence_gauche(floor(essieux_ind)) = 1 - essieux_ind + floor(essieux_ind);
influence_droite = zeros(1, N);
influence_droite(floor(essieux_ind)+1) = essieux_ind - floor(essieux_ind);
influence = influence_gauche + influence_droite;
influence = influence(3:N-1); % supression des ddl des CL
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
