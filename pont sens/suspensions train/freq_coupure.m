% Hunting Stability of High-Speed Railway Vehicles Under Steady Aerodynamic
% Loads, Zeng, 2017
mc = 33786;
k1 = 886500;
k2 = 2.03e5;
K1 = 4*k1;
K2 = 8*k2;
K = K1*K2/(K1+K2);
f = sqrt(K/mc)/(2*pi)

% freq mode vertical

z = 1/sqrt(2);
H = @(x) abs((1 + 2i*z*x)/(1 - x^2 + 2i*z*x));
1/H(4.2/f)