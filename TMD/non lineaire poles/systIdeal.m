%% params
mu = 0.01;
w0 = 2*pi;

x0 = [0 0];
v0 = [1 0];



n = 5;
alphas = linspace(0, 2, n);
lambdas = nan(1, n);
omegas = nan(1, n);
epsilons = nan(1, n);
convergence = false(1, n);


options = optimset('MaxIter',1e4, 'MaxFunEvals', 1e4);

%%

M = @(mu, x0, v0, alpha, X) max (real (InitPolesNonLin(mu, w0, X(1), X(2), alpha, x0, v0)));

w = waitbar(0, '');

X = [1 1];
for ka = 1:n
    alpha = alphas(ka);
    
    waitbar((ka-1)/n, w, ['\alpha = ' num2str(alpha)]);
    
%     [X, fval, exitflag, output] = fmincon(@(x) M(mu, x0, v0, alpha, x), [w0 1], options);
    
    [X, fval, exitflag, output] = fmincon(@(x) M(mu, x0, v0, alpha, x), [w0 0], [1 0;-1 0], w0*[1.5; -0.5]);
    
    lambdas(ka) = M(mu, x0, v0, alpha, X);
    omegas(ka) = X(1);
    epsilons(ka) = X(2);
    
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
plot(ax, alphas(convergence), epsilons(convergence));
plot(ax, alphas(~convergence), epsilons(~convergence), '*r');
hold(ax, 'off');
xlabel(ax, '\alpha');
ylabel(ax, '\epsilon');

