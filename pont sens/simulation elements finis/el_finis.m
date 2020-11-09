clear all

save_results = 0;
results_name = 'muLw1';
results_folder = 'testsMasseAjoutee';

plot_results = 1;

%% données pont

L = 17.5; % longueur du pont
E = 30e9; % module d'young du béton
rho = 2500; % masse volumique du béton
l_pont = 4.85; % largeur du pont
h_pont = 1.2;
J = l_pont * h_pont^3 / 12; % moment quadratique du pont
mu = l_pont * h_pont * rho; % masse linéique du pont
amort = 1.5e4; % coefficient d'amortissement


%% integration spatiale

N = 50; % nb de points ds l'espace
dx = L / (N-1);

%% données capteurs ponts

% pos_capteurs = [L/6, L/3, L/2, 2*L/3, 5*L/6];
% pos_capteurs = linspace(0, L, N); pos_capteurs = pos_capteurs(2:end-1);
pos_capteurs = L/2;


%% données train

for mu_t = [0.004, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]*mu % masse linéique du train
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%');
    clockTime = clock;
    fprintf('mu_t = %.2fmu (%d:%02d)\n', [mu_t/mu, clockTime(4:5)]);
    
    L_wagons = 18.5; % écart entre les wagons
%     L_wagons = 1; % test
    N_bogies = floor(300/L_wagons)+1; % nombre de bogies
    c = 78.5; % vitesse du train
    g = 9.81;
    
    %test
    N_bogies = ceil(200*c/L_wagons) + 1;
    
    %% données temps
    
    ti = 0;
    tf = 50;
    dt = 0.01;
    
    t = ti:dt:tf;
    
    t1 = 0;
    t2 = ((N_bogies-1)*L_wagons + L)/c;
    disp(sprintf('t1 = %.2f, t2 = %.2f', [t1, t2]));
    
    
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
    K_CL1 = 1 * K_CL1;
    
    K_CL2 = zeros(N); % CL \partial^2\nu/\partial x^2(x=0) = \partial^2\nu/\partial x^2(x=L) = 0
    K_CL2(2, 1:3) = [1 -2 1];
    K_CL2(N-1, N-2:N) = [1 -2 1];
    K_CL2 = 1/dx^2 * K_CL2;
    
    K = K_eq + K_CL1 + K_CL2;
    
    % amortissement
%     C = zeros(N);
%     for i = 3:N-2
%         C(i, i-1:i+1) = [1 -2 1];
%     end
%     C = amort/dx^2 * C;
    C = amort * diag([0, 0, ones(1, N-4), 0, 0]);
    
    % opérateur accélération coriolis
    A1 = zeros(N);
    for i = 3:N-2
        A1(i, i-1:i+1) = [-1/2 0 1/2];
    end
    A1 = 2*c/dx * A1;
    
    % opérateur accélération centripète
    A2 = zeros(N);
    for i = 3:N-2
        A2(i, i-1:i+1) = [1 -2 1];
    end
    A2 = c^2/dx^2 * A2;
    
    % restriction aux ddl 3 à N-1 (sinon M non inversible)
    % (y1 = yN = 0, et y2 = y3/2 et yN-1 = yN-2/2)
    Mtransition = diag([0, 0, ones(1, N-4), 0, 0]);
    Mtransition(2, 3) = 1/2;
    Mtransition(N-1, N-2) = 1/2;
    
    % matrices réduites
    Mr = M * Mtransition;
    Cr = C * Mtransition;
    Kr = K * Mtransition;
    A1r = A1 * Mtransition;
    A2r = A2 * Mtransition;
    Mr = Mr(3:N-2, 3:N-2);
    Cr = Cr(3:N-2, 3:N-2);
    Kr = Kr(3:N-2, 3:N-2);
    A1r = A1r(3:N-2, 3:N-2);
    A2r = A2r(3:N-2, 3:N-2);
    
    %% calcul freq propre, freq excitation
    
    % frequences therorique, problème continu
    freqsTh = pi/(2*L^2) * sqrt(E*J/mu) * (1:N-4).^2;
    
    % frequence propre problème discrétisé
    freqs = eig(Mr\Kr);
    freqs = sqrt(freqs)/(2*pi);
    freqs = sort(freqs);
    
    % affichage
    for kfreq = 1:1
        fprintf('freq. propre %d : %.2fHz, th. %.2fHz (%.1f%% error)\n',...
            [kfreq, freqs(kfreq), freqsTh(kfreq), (freqs(kfreq)-freqsTh(kfreq))/freqsTh(kfreq)*100]);
    end
    
    fprintf('freq. excitation : %.2f\n', c/L_wagons);
    
    
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
    
    if any(diff(essieux) <= 2*dx)
        warning('écart d''essieux < 2dx');
    end
    
    % %test
    % r = 0.25;
    % f0 = zeros(N-4, 1);
    % f0(floor(r*(N-1))+1-2) = 1 - r*(N-1) + floor(r*(N-1));
    % f0(floor(r*(N-1))+1-2+1) = r*(N-1) - floor(r*(N-1));
    % F = @(t) - 1e8/dx * f0;
    
    
    %% integration numerique
    D = @(t, YdY) [ YdY(N-3:end);
        DD(t, YdY, dx, essieux, c, mu_t, L_wagons, N, g, Mr, Cr, Kr, A1r, A2r)];
    
    Y0 = zeros(N-4, 1);
    V0 = zeros(N-4, 1);
    YV0 = [Y0; V0];
    
    % temps de calcul
    computationTime = now;
    
    % integration
    options = odeset('OutputFcn', @waitbarOutput);
    [t, YV] = ode45(D, t, YV0, options);
    
    % temps de calcul
    computationTime = now - computationTime;
    computationTime = days(computationTime);
    formatsTime = {'MM:SS', 'HH:MM:SS'};
    computationTime = datestr(computationTime, formatsTime{logical(hms(computationTime))+1});
    
    % mise en forme des matrices
    t = t';
    YV = YV';
    Y = YV(1:N-4, :);
    
    Ytot = getYtot(Y);
    
    
    %% affichage
    
    if plot_results || true
        % animation
        moving_coeff = 1.;
        movingPlot(Ytot, t, L, essieux, c, pos_capteurs, moving_coeff);
        
        % wavelet capteurs
        Ycapt = getYcapt(Ytot, pos_capteurs, dx);
        %Ycapt = Ycapt + 1e-6 * randn(size(Ycapt));
        
        
        fig = figure;
        ax = axes(fig);
        plt = plot(ax, t, Ycapt);
        
        fmin = 3;
        fmax = 16;
        Q = 10;
        MaxRidges = 1;
        RealShapePlot = deformeePont(L, pos_capteurs);
        WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', false,...
            'MaxRidges', MaxRidges, 'RealShapePlot', RealShapePlot);
    end
    
    
    %% enregistrement
    
    if save_results
        folder_dir = 'pont sens/simulation elements finis/resultats';
        if ~isempty(results_folder)
            folder_dir = [folder_dir, '/', results_folder];
            [~,~,~] = mkdir(folder_dir);
        end
        
        listing = dir(folder_dir);
        listingNames = {listing.name};
        nbSimul = 0;
        for kname = 1:length(listingNames)
            name1 = strsplit(listingNames{kname}, '_');
            name1 = name1{1};
            if length(name1) > 5 && strcmp(name1(1:5), 'simul')
                nbSimul = max(nbSimul, str2double(name1(6:end)));
            end
        end
        if ~isempty(results_name)
            simul_name = [sprintf('simul%d_', nbSimul+1), results_name, sprintf('_N%d', N)];
        else
            simul_name = sprintf('simul%d_N%d', [nbSimul+1, N]);
        end
        save([folder_dir, '/', simul_name],...
            'L', 'E', 'rho', 'l_pont', 'h_pont', 'J', 'mu', 'amort', 'N', 'dx', 'pos_capteurs', 'mu_t',...
            'L_wagons', 'N_bogies', 'c', 'g', 'ti', 'tf', 'dt', 't', 'Ytot', 'Y0', 'V0',...
            'essieux', 'computationTime', 'freqs', 'freqsTh', 't1', 't2');
    end
    
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

%% fonction integratio numerique

function ddY = DD(t, YdY, dx, essieux, c, mu_t, L_wagons, N, g, Mr, Cr, Kr, A1r, A2r)
% influence essieux
essieux_ind = (essieux + c*t) / dx + 1;
essieux_ind = essieux_ind(essieux_ind >= 1 & essieux_ind < N);
influence_gauche = zeros(1, N);
influence_gauche(floor(essieux_ind)) = 1 - essieux_ind + floor(essieux_ind);
influence_droite = zeros(1, N);
influence_droite(floor(essieux_ind)+1) = essieux_ind - floor(essieux_ind);
influence = influence_gauche + influence_droite;
influence = influence(3:N-2); % supression des ddl des CL

% matrice de masse ajoutee
M_ajout = mu_t*L_wagons/dx * diag(influence);

% excitation pesanteur
F = -g * mu_t*L_wagons/dx  * transpose(influence);

% acceleration
ddY = (Mr + M_ajout) \ (- (Kr + M_ajout*A2r) * YdY(1:N-4) - (Cr + M_ajout*A1r) * YdY(N-3:end) + F);
end


%% interpolation etc

function Ytot = getYtot(Y) % rajoute les DDL 1, 2, N-1 et N
Ytot = [zeros(1, size(Y, 2)); Y(1, :)/2; Y; Y(end, :)/2; zeros(1, size(Y, 2))];
end


%% barre de chargement

function status = waitbarOutput(t, ~, flag)
persistent ti tf f t_last t_1s T_1s t0

update_waitbar_time = 0.1;
update_remainingTime_time = 2;
minimum_beep_time = 1*60;

status = 0;
if nargin < 3 || isempty(flag)
    t_0s = 24*3600*now;
    if t_0s - t_last > update_waitbar_time
        t_last = t_0s;
        x = (t(end)-ti)/(tf-ti);
        try
            waitbar(x, f, sprintf('t = %.2f', t(end)));
        catch
            status = 1;
        end
    end
    if t_0s - t_1s > update_remainingTime_time
        t_remaining = (t_0s - t_1s) * (tf-t(end))/(t(end)-T_1s);
        try
            set(f, 'Name', ['Computing ODE (', timeString(t_remaining),' remaining)']);
        catch
            status = 1;
        end
        t_1s = t_0s;
        T_1s = t(end);
    end
elseif strcmp(flag, 'init')
    f = waitbar(0, sprintf('t = %.2f', t(1)), 'Name', 'Computing ODE');
    ti = t(1);
    tf = t(end);
    t_last = 24*3600*now;
    t_1s = t_last;
    t0 = t_last;
    T_1s = t(1);
elseif strcmp(flag, 'done')
    if 24*3600*now - t0 >= minimum_beep_time
        beep
    end
    try
        close(f);
    catch
    end
end
end

function timeS = timeString(t)
[h, m, s] = hms(seconds(t));
s = 10*floor(s/10);
if h == 0 && m == 0
    timeS = sprintf('%ds', s);
elseif h == 0 && m < 5 && s ~= 0
    timeS = sprintf('%dm%02ds', [m, s]);
elseif h == 0
    timeS = sprintf('%dm', m);
else
    timeS = sprintf('%dh%02dm', [h, m]);
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
