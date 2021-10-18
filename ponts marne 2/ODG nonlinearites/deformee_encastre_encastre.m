betaL = 4.73004074;
sig = 0.9825;

L = 65.74;
beta = betaL/L;

phi = @(x) cosh(beta*x) - cos(beta*x) - sig*(sinh(beta*x) - sin(beta*x));
phi2 = @(x) beta^2*(cosh(beta*x) + cos(beta*x) - sig*(sinh(beta*x) + sin(beta*x)));

% l = linspace(0, L);
% figure;
% yyaxis left
% plot(l, phi(l));
% yyaxis right
% plot(l, phi2(l));

phi_mitravee = phi(L/2);
phi2max = phi2(0);

%%

E = 40e9;
sigma0 = 8e6;
e = 2.4;
% e = 0.86;

ampl_ouverture_fissures = sigma0*phi_mitravee/(phi2max*E*e/2)

%%

A_acc_max = 0.29;
f1 = 2.14;
A_max = A_acc_max/(2*pi*f1)^2