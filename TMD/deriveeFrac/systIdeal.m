mu = 0.01;


n = 80;

n = n-1;
alphas = linspace(0, 2, n+2);
alphas = alphas(2:end-1);
n = 1;
alphas = 1;
lambdas = nan(1, n);
omegas = nan(1, n);
zetas = nan(1, n);
convergence = false(1, n);


options = optimset('MaxIter',1e4, 'MaxFunEvals', 1e4);

%%

M = @(mu, alpha, X) max (real (polesSystFrac(mu, 1, X(1), 0, X(2), alpha)));

w = waitbar(0, '');

X = [1 1];
for ka = 1:n
    alpha = alphas(ka);
    
    waitbar((ka-1)/n, w, ['\alpha = ' num2str(alpha)]);
    
    [X, fval, exitflag, output] = fminsearch(@(x) M(mu, alpha, x), [1 1], options);
    
    lambdas(ka) = M(mu, alpha, X);
    omegas(ka) = X(1);
    zetas(ka) = X(2);
    
    convergence(ka) = exitflag;
    
    waitbar(ka/n, w);
end
close(w);

%%

f = figure;
ax = axes(f);
hold(ax, 'on');
plot(ax, alphas(convergence), lambdas(convergence));
plot(ax, alphas(~convergence), lambdas(~convergence), '*r');
hold(ax, 'off');
xlabel(ax, '\alpha');
ylabel(ax, '-\lambda');


f = figure;
ax = axes(f);
hold(ax, 'on');
plot(ax, alphas(convergence), omegas(convergence));
plot(ax, alphas(~convergence), omegas(~convergence), '*r');
hold(ax, 'off');
xlabel(ax, '\alpha');
ylabel(ax, '\omega');


f = figure;
ax = axes(f);
hold(ax, 'on');
plot(ax, alphas(convergence), zetas(convergence));
plot(ax, alphas(~convergence), zetas(~convergence), '*r');
hold(ax, 'off');
xlabel(ax, '\alpha');
ylabel(ax, '\zeta');



RegressionMenu;
