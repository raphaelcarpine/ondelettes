beta_L_sur_pi = 0.596864;

beta_L = beta_L_sur_pi*pi;

coeff_freq = beta_L^2/(2*pi)


L = 1;
x = linspace(0, L, 100000);
w = cosh(beta_L*x/L) - cos(beta_L*x/L) + (cos(beta_L)+cosh(beta_L))/(sin(beta_L)+sinh(beta_L)) * (sin(beta_L*x/L) - sinh(beta_L*x/L));
w = w/w(end);

coeff = mean(w).^2/mean(w.^2)