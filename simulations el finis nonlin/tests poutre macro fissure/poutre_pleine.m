% https://linkinghub.elsevier.com/retrieve/pii/S0022460X98916406

zeta = 0.3/2.4; % a/b, avec a profondeur crack et b épaisseur totale

phi1 = @(z) 0.6272*z.^2 - 1.04533*z.^3 + 4.5948*z.^4 - 9.9736*z.^5 + 20.2948*z.^6 - 33.0351*z.^7 + 47.1063*z.^8 - 40.7556*z.^9 + 19.6*z.^10;

Kr = E*J/(6*pi*(1-0.2^3)*h_pont*phi1(zeta));

%%

Kr0 = E*J/dx;

Kr2 = 1/(1/Kr + 1/Kr0);

slope = Kr2/Kr0