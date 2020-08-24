%% données pont

L = 17.5; % longueur du pont
E = 30e9; % module d'young du béton
rho = 2500; % masse volumique du béton
l_pont = 4.85; % largeur du pont
h_pont = 1;
J = l_pont * h_pont^3 / 12; % moment quadratique du pont
mu = l_pont * h_pont * rho; % masse linéique du pont
amort = 4e5; % coefficient d'amortissement

%% capteurs ponts

pos_capteurs = [L/6, L/3, L/2, 2*L/3, 5*L/6];


%% données train

mu_t = mu * 0.5; % masse linéique du train
L_bogies = 18.5; % écartement entre les bogies
L_bogies = 3; % écartement entre les bogies
l_bogies = 2; % écartement des deux essieux au sein d'un bogie
l_bogies = 10;
N_wagons = 16; % nombre de wagons
N_wagons = 1;
c = 78.5; % vitesse du train
g = 9.81;


%% données temps

ti = -2;
tf = 10;
% t0 = 0;
dt = 0.01;


%% integration spatiale

N = 80; % nb de points ds l'espace
dx = L / (N-1);
%x = linspace(0, L, N);

if any([l_bogies, l_bogies, l_bogies-l_bogies] <= 2*dx)
    warning('');
end


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
M_ajout = @(t) mu_t * L_bogies/2 * diag( influence_train(essieux + c*t, dx, N));

% excitation pesanteur
F = @(t) -g*mu_t * L_bogies/2  * transpose( influence_train(essieux + c*t, dx, N));


%% integration numerique

D = @(t, YdY) [ YdY(N-3:end);
    (Mr + M_ajout(t)) \ (-Kr * YdY(1:N-4) - Cr * YdY(N-3:end) + F(t))];

X0 = zeros(N-4, 1);
V0 = zeros(N-4, 1);
XV0 = [X0; V0];

[t, XV] = ode45(D, [ti, tf], XV0);
t = t';
XV = XV';
X = XV(1:N-4, :);

Xtot = getXtot(X);


%% interpolation

t_interp = ti:dt:tf;
Xtot_interp = interp1(t, Xtot', t_interp)';

%% affichage

% animation
movingPlot(Xtot_interp, t_interp, L, essieux, c, pos_capteurs);

% wavelet capteurs
Xcapt_interp = getXcapt(Xtot_interp, pos_capteurs, N, dx);
Xcapt_interp = Xcapt_interp + 1e-6 * randn(size(Xcapt_interp));


fig = figure;
ax = axes(fig);
plt = plot(ax, t_interp, Xcapt_interp);

fmin = 3;
fmax = 16;
Q = 10;
MaxRidges = 6;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', true, 'MaxRidges', MaxRidges);




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



%%

function Xtot = getXtot(X)
Xtot = [zeros(1, size(X, 2)); X(1, :)/2; X; X(end, :)/2; zeros(1, size(X, 2))];
end

function Xcapt = getXcapt(Xtot, pos_capt, N, dx)
Xcapt = nan(length(pos_capt), size(Xtot, 2));
for it = 1:size(Xtot, 2)
    Xcapt(:, it) = transpose(...
        Xtot( max(min( floor(pos_capt/dx) + 1, N), 1), it)' .* (1 - pos_capt/dx + floor(pos_capt/dx))...
        + Xtot( max(min( floor(pos_capt/dx) + 2, N), 1), it)' .* (pos_capt/dx - floor(pos_capt/dx)));
end
end


%% animation

function movingPlot(Xtot, t, L, essieux, c, pos_capteurs)
timeCoeff = 1.;

% framerate interpolation
framerate = 20;
t_frames = t(1):(timeCoeff/framerate):t(end);
Xtot = interp1(t, Xtot', t_frames)';
t = t_frames;

N = size(Xtot, 1);
dx = L / (N-1);
Xtot = Xtot / max(abs(Xtot), [], 'all');

fig = figure;
ax = axes(fig);
hold(ax, 'on');

xlim(ax, [-L/3, 4*L/3]);
ylim(ax, [-5, 5]);
xticks([0, L]);
xticklabels({'0', 'L'});
yticks([]);

    function animation()
        cla(ax);
        
        plt1 = plot(ax, linspace(0, L, N), Xtot(:, 1), '-+');
        essieux_1 = essieux + c*t(1);
        plt2 = scatter(ax, essieux_1,...
            Xtot( max(min( floor(essieux_1/dx) + 1, N), 1), 1)' .* (1 - essieux_1/dx + floor(essieux_1/dx))...
                + Xtot( max(min( floor(essieux_1/dx) + 2, N), 1), 1)' .* (essieux_1/dx - floor(essieux_1/dx)), 'r');
        plt3 = scatter(ax, pos_capteurs,...
            Xtot( max(min( floor(pos_capteurs/dx) + 1, N), 1), 1)' .* (1 - pos_capteurs/dx + floor(pos_capteurs/dx))...
                + Xtot( max(min( floor(pos_capteurs/dx) + 2, N), 1), 1)' .* (pos_capteurs/dx - floor(pos_capteurs/dx)), 'dm');
        txt = text(ax, 0, 4, ['t = ', num2str(t(1))]);
        
        pause(1);
        
        tic;
        for it = 1:length(t)
            set(plt1, 'YData', Xtot(:, it));
            essieux_it = essieux + c*t(it);
            set(plt2, 'XData', essieux_it, 'YData',...
                Xtot( max(min( floor(essieux_it/dx) + 1, N), 1), it)' .* (1 - essieux_it/dx + floor(essieux_it/dx))...
                + Xtot( max(min( floor(essieux_it/dx) + 2, N), 1), it)' .* (essieux_it/dx - floor(essieux_it/dx)));
            set(plt3, 'YData',...
                Xtot( max(min( floor(pos_capteurs/dx) + 1, N), 1), it)' .* (1 - pos_capteurs/dx + floor(pos_capteurs/dx))...
                + Xtot( max(min( floor(pos_capteurs/dx) + 2, N), 1), it)' .* (pos_capteurs/dx - floor(pos_capteurs/dx)));
            set(txt, 'String', ['t = ', num2str(t(it))]);
            drawnow
            if it < length(t)
                pause( t(it+1) - t(1) - toc*timeCoeff);
            end
        end
    end

%animation();

set(ax, 'ButtonDownFcn', @(~, ~) animation);

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
