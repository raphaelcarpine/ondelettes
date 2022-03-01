% https://doi.org/10.1016/0167-8442(90)90023-S

zeta = 1/6; % a/b, avec a profondeur crack et b épaisseur totale

Yf = @(z) 1.99*z.^(1/2) - 0.41*z.^(3/2) + 18.70*z.^(5/2) - 38.48*z.^(7/2) + 53.85*z.^(9/2);
z0 = linspace(0, zeta, 10000);
dz = mean(diff(z0));
pz = Yf(z0).^2;
pz = dz/2 * sum(pz(1:end-1) + pz(2:end));

% K = E*l_pont/(2*pz); % fissure pas sur toute l'épaisseur, mais toute la largeur
K = E*e_pont/(2*pz); % fissure sur toute l'épaisseur, mais pas toute la largeur

Kr = (h_pont-e_pont)^2*K;

%%

Kr0 = E*J/dx;

Kr2 = 1/(1/Kr + 1/Kr0);

slope = Kr2/Kr0