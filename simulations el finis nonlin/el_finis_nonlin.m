function el_finis_nonlin()

results_name = 'test';
saveFolder = 'simulations el finis nonlin\resultats';

solve_ODE = 1;
plot_results = 1;

nonLin = 1;
inertia_vehicles = 0;
shock_mode = 1;

%% données pont

L = 70; % longueur du pont
E = 30e9; % module d'young du béton
rho = 2500; % masse volumique du béton
l_pont = 8; % largeur du pont
h_pont = 3;
J = l_pont * h_pont^3 / 12; % moment quadratique du pont
mu = l_pont * h_pont * rho; % masse linéique du pont
amort_deg = 2; % degré derivation spatial amortissement
amort = -3.75e6; % coefficient d'amortissement, 7.6e3, -3.75e6, 1.71e9


%% integration spatiale

N = 50; % nb de points ds l'espace
dx = L / (N-1);

%% non linearite

x_nonlin = L/pi;
px_nonlin = round(x_nonlin/dx)+1;
x_nonlin = (px_nonlin-1)*dx;
threshold_nonlin = 0; % courbure
slope_nonlin = 0.1; % *K

%% données capteurs ponts

% pos_capteurs = [L/6, L/3, L/2, 2*L/3, 5*L/6];
pos_capteurs = linspace(0, L, 21);% pos_capteurs = pos_capteurs(2:end-1);
% pos_capteurs = L/2;


%% données temps

t_tot = 100;
ti = 0;
tf = t_tot;
Fe = 2000;
dt = 1/Fe;

T = ti:dt:tf;


%% données véhicules

disp('%%%%%%%%%%%%%%%%%%%%%%%%%');
clockTime = clock;
fprintf('%d:%02d\n', clockTime(4:5));

g = 9.81;

tau_vehicle = [1 1]; % constante de tps loi sans memoire, deux sens de circulation
m0_vehicle = 1e3; % espérance masse véhicules
sigma_m_vehicle = 3e2; % écart type masse véhicules
c0_vehicle = 60/3.6; % espérance vitesse véhicules
sigma_c_vehicle = 10/3.6; % écart type vitesse véhicules

if shock_mode
    tau_vehicle = inf*tau_vehicle; % pas de vehicule
end

% tirages aléatoire trafic
t_vehicles_left = 0;
while t_vehicles_left(end) <= t_tot
    t_vehicles_left(end+1) = t_vehicles_left(end) + exprnd(tau_vehicle(1));
end
t_vehicles_left = t_vehicles_left(2:end-1);
m_vehicles_left = m0_vehicle + sigma_m_vehicle*randn(size(t_vehicles_left));
c_vehicles_left = c0_vehicle + sigma_c_vehicle*randn(size(t_vehicles_left));

t_vehicles_right = 0;
while t_vehicles_right(end) <= t_tot
    t_vehicles_right(end+1) = t_vehicles_right(end) + exprnd(tau_vehicle(end));
end
t_vehicles_right = t_vehicles_right(2:end-1);
m_vehicles_right = m0_vehicle + sigma_m_vehicle*randn(size(t_vehicles_right));
c_vehicles_right = c0_vehicle + sigma_c_vehicle*randn(size(t_vehicles_right));


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
if amort_deg == 0
    C = amort * diag([0, 0, ones(1, N-4), 0, 0]);
elseif amort_deg == 2
    C = zeros(N);
    for i = 3:(N-2)
        C(i, i-1:i+1) = [1 -2 1];
    end
    C = amort/dx^2 * C;
elseif amort_deg == 4
    C = zeros(N);
    for i = 3:(N-2)
        C(i, i-2:i+2) = [1 -4 6 -4 1];
    end
    C = amort/dx^4 * C;
end

% opérateur accélération coriolis
A1 = zeros(N);
for i = 3:N-2
    A1(i, i-1:i+1) = [-1/2 0 1/2];
end
A1 = 2/dx * A1; % à multiplier par la vitesse

% opérateur accélération centripète
A2 = zeros(N);
for i = 3:N-2
    A2(i, i-1:i+1) = [1 -2 1];
end
A2 = 1/dx^2 * A2; % à multiplier par la vitesse au carré

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

% precalcul pour cas lineaire dans inertie vehocules
Mr_Kr = Mr\Kr;
Mr_Cr = Mr\Cr;
invMr = inv(Mr);

% derivation nonlin
d2_nonlin = zeros(1, N-4);
d2_nonlin(px_nonlin-3:px_nonlin-1) = 1/dx^2 * [1 -2 1];
D2_nonlin = A2r; % derivation^2
courbure_vect0 = zeros(N-4, 1);
courbure_vect0(px_nonlin-2) = 1;

%% calcul freq propre, freq excitation

% frequences therorique, problème continu
freqsTh = pi/(2*L^2) * sqrt(E*J/mu) * (1:N-4).^2;

% frequence propre problème discrétisé
freqs = eig(Mr\Kr);
freqs = sqrt(freqs)/(2*pi);
freqs = sort(freqs);

% poles
poles = eig([zeros(size(Mr)), eye(size(Mr)); -Mr\Kr, -Mr\Cr]);
poles = poles(imag(poles) ~= 0);
[~, Ipoles] = sort(imag(poles));
poles = poles(Ipoles(end/2+1:end));
amorts = -real(poles)./abs(poles);
amorts = [amorts.', inf*ones(1, length(freqs)-length(amorts))];

% affichage
for kfreq = 1:5
    fprintf('mode %d: f=%.2fHz (f_th=%.2fHz, %.1f%% error), z=%.2f%%\n',...
        [kfreq, freqs(kfreq), freqsTh(kfreq), (freqs(kfreq)-freqsTh(kfreq))/freqsTh(kfreq)*100, 100*amorts(kfreq)]);
end
disp('...');
kfreq = length(freqs);
fprintf('mode %d: f=%.2fHz (f_th=%.2fHz, %.1f%% error), z=%.2f%%\n',...
    [kfreq, freqs(kfreq), freqsTh(kfreq), (freqs(kfreq)-freqsTh(kfreq))/freqsTh(kfreq)*100, 100*amorts(kfreq)]);


if ~solve_ODE
    return
end

%% integration numerique
D = @(t, YdY) [ YdY(N-3:end);
    DD(t, YdY)];

Y0 = zeros(N-4, 1);
V0 = zeros(N-4, 1);
if shock_mode
    V0 = rand(N-4, 1); % choc
end
YV0 = [Y0; V0];

% temps de calcul
computationTime = now;

% wait bar
[initWaitBar, updateWaitBar, closeWaitBar] = getWaitBar(length(T), 'windowTitle', 'Computing ODE',...
    'msg', '', 'minimumBeepTime', 0);
initWaitBar();
    function flag = updateWaitBar2(t, ~, ~)
        if ~isempty(t)
            updateWaitBar((t(1)-T(1))/mean(diff(T)), sprintf('t = %.1fs', t(1)))
        end
        flag = 0;
    end

% integration
options = odeset('OutputFcn', @updateWaitBar2);
[T, YV] = ode45(D, T, YV0, options);
closeWaitBar();

% temps de calcul
computationTime = now - computationTime;
computationTime = days(computationTime);
formatsTime = {'MM:SS', 'HH:MM:SS'};
computationTime = datestr(computationTime, formatsTime{logical(hms(computationTime))+1});

% mise en forme des matrices
T = T';
YV = YV';
Y = YV(1:N-4, :);
V = YV(N-3:end, :);

% moments où la nonlinearite est atteinte
nonlin_reached = abs(d2_nonlin * Y) >= threshold_nonlin;
nonlin_reached = nonlin_reached & nonLin;

A = nan(size(YV));
for kt = 1:length(T)
    A(:, kt) = D(T(kt), YV(:, kt));
end
A = A(N-3:end, :);

Ytot = getYtot(Y);
Vtot = getYtot(V);
Atot = getYtot(A);


%% enregistrement

listing = dir(saveFolder);
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
    simul_name = [sprintf('simul%d_', nbSimul+1), results_name, sprintf('_N%d_T%d', [N, t_tot])];
else
    simul_name = sprintf('simul%d_N%d_T%d', [nbSimul+1, N, t_tot]);
end
%     save([folder_dir, '/', simul_name],...
%         'L', 'E', 'rho', 'l_pont', 'h_pont', 'J', 'mu', 'amort', 'N', 'dx', 'pos_capteurs', 'mu_t',...
%         'L_wagons', 'N_bogies', 'c', 'g', 'ti', 'tf', 'dt', 't', 'Ytot', 'Y0', 'V0',...
%         'essieux', 'computationTime', 'freqs', 'freqsTh', 't1', 't2');
save(fullfile(saveFolder, simul_name));

disp('saved');


%% affichage

if plot_results
    el_finis_nonlin_display();
end


%% fonction integratio numerique


    function ddY = DD(t, YdY)
        % vehicules sur le pont
        I_left = t_vehicles_left <= t & t_vehicles_left + L./c_vehicles_left >= t;
        I_right = t_vehicles_right <= t & t_vehicles_right + L./c_vehicles_right >= t;
        pos_left = c_vehicles_left(I_left) .* (t - t_vehicles_left(I_left));
        pos_right = L - c_vehicles_right(I_right) .* (t - t_vehicles_right(I_right));
        pos_on_bridge = [pos_left, pos_right];
        m_on_bridge = [m_vehicles_left(I_left), m_vehicles_right(I_right)];
        
        % influence poids
        pos_ind = pos_on_bridge / dx + 1;
        F = zeros(N, 1);
        for k_v = 1:length(pos_ind)
            % séparation sur les deux points les plus proches
            p0 = floor(pos_ind(k_v));
            F(p0) = F(p0) + (p0+1 - pos_ind(k_v)) * m_on_bridge(k_v); % bras de levier
            F(p0+1) = F(p0+1) + (pos_ind(k_v) - p0) * m_on_bridge(k_v);
        end
        
        F = - g/dx*F(3:N-2);
        
        if inertia_vehicles
            c_on_bridge = [c_vehicles_left(I_left), -c_vehicles_right(I_right)];
            
            % influence  inertie
            pos_ind = pos_on_bridge / dx + 1;
            pos_ind = round(pos_ind); % pas de séparation sur les deux points pour l'inertie
            influenceM = zeros(1, N);
            influenceMc = zeros(1, N);
            influenceMc2 = zeros(1, N);
            for k_v = 1:length(pos_ind)
                influenceM(pos_ind(k_v)) = influenceM(pos_ind(k_v)) + m_on_bridge(k_v);
                influenceMc(pos_ind(k_v)) = influenceMc(pos_ind(k_v)) + m_on_bridge(k_v) * c_on_bridge(k_v);
                influenceMc2(pos_ind(k_v)) = influenceMc2(pos_ind(k_v)) + m_on_bridge(k_v) * c_on_bridge(k_v)^2;
            end
            
            % matrice de masse ajoutee
            M_ajout = 1/dx * diag(influenceM(3:N-2)); % supression des ddl des CL
            M_ajout_c = 1/dx * diag(influenceMc(3:N-2));
            M_ajout_c2 = 1/dx * diag(influenceMc2(3:N-2));
        elseif ~inertia_vehicles && nonLin
            % matrice de masse ajoutee
            M_ajout = zeros(N-4);
            M_ajout_c = zeros(N-4);
            M_ajout_c2 = zeros(N-4);
        end
        
        if ~inertia_vehicles && ~nonLin
            % acceleration
            ddY = -Mr_Kr * YdY(1:N-4) - Mr_Cr * YdY(N-3:end) + invMr * F;
        elseif inertia_vehicles && ~nonLin
            ddY = (Mr + M_ajout) \ (- (Kr + M_ajout_c2*A2r) * YdY(1:N-4) - (Cr + M_ajout_c*A1r) * YdY(N-3:end) + F);
        else % nonLin
            courbure = d2_nonlin * YdY(1:N-4);
            if courbure > threshold_nonlin
                F_nonlin = (1-slope_nonlin)*E*J * D2_nonlin * (courbure-threshold_nonlin)*courbure_vect0;
            elseif courbure < -threshold_nonlin
                F_nonlin = (1-slope_nonlin)*E*J * D2_nonlin * (courbure+threshold_nonlin)*courbure_vect0;
            else
                F_nonlin = zeros(N-4, 1);
            end
            
            ddY = (Mr + M_ajout) \ (- (Kr + M_ajout_c2*A2r) * YdY(1:N-4) - (Cr + M_ajout_c*A1r) * YdY(N-3:end) + F + F_nonlin);
        end
        
        
    end

%% interpolation etc

    function Ytot = getYtot(Y) % rajoute les DDL 1, 2, N-1 et N
        Ytot = [zeros(1, size(Y, 2)); Y(1, :)/2; Y; Y(end, :)/2; zeros(1, size(Y, 2))];
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

end