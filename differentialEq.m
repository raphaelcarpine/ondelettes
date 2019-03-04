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

switch visible
    case true
        stats = 'on';
    case false
        stats = 'off';
end

[t, X] = solver(f, [0 T], x0, odeset('RelTol',1e-10,'Stats',stats)); %,'OutputFcn',@odeplot
n = size(X,2)/2;
x = X(:, 1:n);
v = X(:, n+1:2*n);
va = zeros(size(X));
for it = 1:length(t)
    va(it, :) = f(t(it), X(it, :).').';
end
a = va(:, n+1:2*n);

if visible
    multiplePlot(t, x, "x");
    multiplePlot(t, v, "v");
    multiplePlot(t, a, "a");
end

    function multiplePlot(t, x, xname)
        n = size(x, 2);
        fig = figure;
        points = {};
        for k=1:n
            subplot(n, 1, k);
            set(gca, 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
            hold on
            plot(t, x(:,k), 'Color', [0 0.4470 0.7410]);
            points{end+1} = plot(t, x(:,k), 'o', 'Color', [0 0.4470 0.7410], 'Visible', 'off');
            hold off
            xlabel('t');
            ylabel(xname + k);
        end
        set(fig, 'WindowButtonDownFcn', @(~, ~) tooglePoints(points));
    end

    function tooglePoints(points)
        toggle = isequal(get(points{1}, 'Visible'), 'off');
        for p = [points{:}]
            set(p, 'Visible', toggle);
        end
    end

end