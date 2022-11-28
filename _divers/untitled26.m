P = 15000; % pret
m = 279.27; % mensualité
T = 60; % temps (mois)

%% cout

C = m*T - P

%% taux

alpha0 = 0.05/12; % taux

F = @(a) m/a + (P-m/a)*exp(a*T);
alpha = 12*lsqnonlin(F, alpha0)