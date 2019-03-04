clear
close all
%% Systeme

%
%   /_/_/_/    ___     ___     ___     ___
%   |   |       |       |       |       |
%   k1  c1      |       |       |       |
%   |   |       | x1    |       |       |
%     T         |       |       |       |
%     m1     _ _V_      |       |       |
%     |                 |       |       |
%   _____               | x2    |       |
%   |   |               |       |       |
%   k2  c2              |       |       |
%   |   |               |       |       |
%     T                 |       | x3    |
%     m2     _ _ _ _ _ _V_      |       |
%     |                         |       |
%   _____                       |       |
%   |   |                       |       | x4
%   k3  c3                      |       |
%   |   |                       |       |
%     T                         |       |
%     m3     _ _ _ _ _ _ _ _ _ _V_      |
%     |                                 |
%   _____                               |
%   |   |                               |
%   k4  c4                              |
%   |   |                               |
%     T                                 |
%     m4     _ _ _ _ _ _ _ _ _ _ _ _ _ _V_
%
%
%
% 8 ddl : 4ddl + 4ddl supp. pour rester avec une EDP du 1er ordre tq
% x"+ x'+ x = f <=> x'= z ; z'+ z + x = f

% On a donc mi zi' = -ki (xi- xi-1) + ki+1 (xi+1 - xi) - ci (xi- xi-1) + ci+1 (xi+1 - xi)
% d/dt xi = zi

%           mi zi' = -ki (xi - xi-1) + ki+1 (xi+1 - xi) - ci zi + ci+1 zi+1
% on choisit [x1;z1;x2;z2;x3;z3;x4;z4]
%
%% Paramètres du système

k1 = 7*1e3;
k2 = 8*1e3;
k3 = 7*1e3;
k4 = 8*1e3;
m1 = 1;
m2 = 1;
m3 = 1;
m4 = 1;
c1 = .7;
c2 = .8;
c3 = .7;c3=1e1;
c4 = .8;

KSys = [k1,m1,c1;k2,m2,c2;k3,m3,c3;k4,m4,c4];

Fun = SysLin4ddlEquaFun(KSys); %
%% Résolution de l'ODE
tf = 50; % Instant final de calcul, on démarre à 0.

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[t0,y0]=ode45(Fun,[0,tf],[0;0;0;0;0;0;0;1],options);
% figure
% plot(t0,y0(:,1:2:end))
%% Interpolation
%
Fs = 200; % Fréquence d'échantillonage souhaitée

t=0:1/Fs:tf;
y = interp1(t0,y0,t);
clear t0 y0
%% Tracé CWT
WvltPlot(t,y(:,7),linspace(1,30,200),25,'PlotScale','log','Title','x_4')
%% Résolution analytique du système : fréquences et amortissements des modes propres
[phi,freq,xi]=SysLin4ddlSol(KSys);
phi = flip(phi,2);
freq = flip(freq);
xi = flip(xi);

%% Identification de ridges
%% Paramètres d'identification
Q = [8,20,25,25]; % Liste des Q retenus par mode
f_min = [4,12,21,25.5] ; % Freq. min d'identification retenues par mode
f_max = [5,15,22,26.5] ; % Freq. max d'identification retenues par mode

MinModu = 1e-8; %Module minimum pour un ridge

%% Paramètres des tracés
PlotDiffExactSup = .25; % lim. sup de la fenêtre de tracé de la grandeur identifiée autour de la valeur exacte, en ratio de la valeur exacte
PlotDiffExactInf = 1; % Lim. inf de la fenêtre de tracé de la grandeur identifiée autour de la valeur exacte,, en ratio de la valeur exacte

% Ex : pour PlotDiffExactSup = .25 et PlotDiffExactInf = .5;
% on trace la freq. identifiée sur un axe allant de f0*(1-.5) à f0*(1+.25)
% f0 étant la valeur théorique.

AmpFactor = .25; % Facteur d'amplification pour le tracé du mode en .gif = valeur de l'amplitude max.
%% Creation (si besoin) du dossier /pics/ à côté de ce .m
p=mfilename('fullpath');
p=p(1:end-length(mfilename));
p=[p,'pics'];
[~,~,~] = mkdir(p);
%% On lance l'extraction et les tracés
RidgeExtractSysLin4ddl