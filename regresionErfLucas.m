%% construction des "données expérimentales"

x = transpose(linspace(0, 10, 100)); % vecteur position (colone)
t = linspace(10, 100, 100); % vecteur temps (ligne)

C0 = 1;
Cdelta = 0.1;
DeltaX = 3;
D0 = 1; % D0 = ke*kt*Drcm(t0)
t0 = 11;
Alpha = 2;

C = C0 + (Cdelta + C0) * (1 - erf((x - DeltaX) ./ (2* sqrt(D0 * (t0./t).^Alpha .* t)))); % forme théorique
C = C + 0.1*randn(size(C)); % bruit ajouté

%% affichage de la surface

figure;
[X, T] = meshgrid(x, t);
surf(X, T, C, 'FaceColor', 'blue', 'FaceAlpha', 0.5);
xlabel('x');
ylabel('t');
zlabel('c');


%% régression non linéaire

% fonction fit
Creg = @(p) p(1) + (p(2) + p(1)) * (1 - erf((x - p(3)) ./ (2* sqrt(p(4) * (t0./t).^p(5) .* t))));
% p = [C0, Cdelta, DeltaX, D0, alpha]

% paramètres initiaux de la régression
C0i = 2;
Cdeltai = 0.01;
DeltaXi = 5;
D0i = 2;
Alphai = 1.5;
Params0 = [C0i, Cdeltai, DeltaXi, D0i, Alphai];

% regression non linéaire
S = @(Params) sum((Creg(Params) - C).^2, 'all'); % somme des carrésà minimiser
ParamsReg = fminsearch(S, Params0); % minimisation


% afffichage
hold on
surf(X, T, Creg(ParamsReg), 'FaceColor', 'red');


c0 = ParamsReg(1)
cdelta = ParamsReg(2)
deltaX = ParamsReg(3)
d0 = ParamsReg(4)
alpha = ParamsReg(5)

