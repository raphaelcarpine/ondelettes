
% fonction de différence entre les poles trouvés avec la transformée en
% ondelettes et les poles théoriques pour certains paramètres
S = @(param) abs(poles - getPoles3ddl(param(1), param(2), param(3), param(4), param(5), param(6)));


f10 = nan; % à compléter
f20 = nan; % à compléter
mu10 = nan; % à compléter
mu20 = nan; % à compléter
ft0 = nan; % à compléter
zetat0 = nan; % à compléter


param0 = [f10, f20, mu10, mu20, ft0, zetat0];


param0 = [6.7, 5.7, 0.015, 0.015, 6.7, 0.07];




% régression permettant de trouver les paramètres minimisant la différence
% entre poles expérimentaux et théoriques
optionsReg = optimoptions(@lsqnonlin, 'MaxIterations', 1e5,...
    'StepTolerance', 1e-6, 'MaxFunctionEvaluations', inf, 'FunctionTolerance', 0);
lbound = zeros(1, 6);
ubound = inf*ones(1, 6);

param = lsqnonlin(S, param0, lbound, ubound, optionsReg);



f1 = param(1);
f2 = param(2);
mu1 = param(3);
mu2 = param(4);
ft = param(5);
zetat = param(6);


% affichage des résultats
disp(['f1 = ', num2str(f1)]);
disp(['f2 = ', num2str(f2)]);
disp(['ft = ', num2str(ft)]);
disp(['zetat = ', num2str(100*zetat), ' %']);
disp(['mu1 = ', num2str(100*mu1), ' %']);
disp(['mu2 = ', num2str(100*mu2), ' %']);
