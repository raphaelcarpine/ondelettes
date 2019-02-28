function [t, x] = differentialEq(x0, f, T, visible, solver)
%differentialEq integre une eq diff du type x'=f(t, x)
%   donne la réponse temporelle sous forme graphique et vectorielle
%   f : fonction
%   x0 : connditions initiales
%   tspan : [ti, tf]
%   solver : type de solveur, 
defaultSolver = @ode45;
switch nargin
    case 3
        visible = true;
        solver = defaultSolver;
    case 4
        solver = defaultSolver;
end

n = uint16(length(x0)/2);
names = repelem("x", length(x0));
for i=1:n
    names(i) = sprintf("x%d", i);
    names(n+i) = sprintf("v%d", i);
end

switch visible
    case true
        stats = 'on';
    case false
        stats = 'off';
end
[t, x] = solver(f, [0 T], x0, odeset('RelTol',1e-6,'Stats',stats)); %,'OutputFcn',@odeplot
if visible
    figure;
    n = size(x,2)/2;
    for k=1:n
        subplot(n, 1, k);
        plot(t, x(:,k), 'o-');
        xlabel('t');
        ylabel(names(k));
    end
end
end