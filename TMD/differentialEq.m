function [t, xout] = differentialEq(x0, f, T, visible, varargin)
%differentialEq integre une eq diff du type x'=f(t, x)
%   donne la réponse temporelle sous forme graphique et vectorielle
%   f : fonction
%   x0 : connditions initiales
%   tspan : [ti, tf]
%   solver : type de solveur, 
p = inputParser;

defaultSolver = @ode45;
defaultOutput = 'x';
validOutputs = {'x', 'v', 'a'};
checkOutput = @(x) ismember(x,validOutputs);
defaultRelTol = 1e-10;
defaultMaxStep = 0.01*T;

addRequired(p,'x0');
addRequired(p,'f');
addRequired(p,'T');
addOptional(p,'visible', true);
addParameter(p,'solver', defaultSolver);
addParameter(p,'output', defaultOutput, checkOutput);
addParameter(p,'RelTol', defaultRelTol);
addParameter(p,'MaxStep', defaultMaxStep);


parse(p, x0, f, T, visible, varargin{:})

solver = p.Results.solver;
output = p.Results.output;
RelTol = p.Results.RelTol;
MaxStep = p.Results.MaxStep;

% switch nargin
%     case 3
%         visible = true;
%         solver = defaultSolver;
%     case 4
%         solver = defaultSolver;
% end

switch visible
    case true
        stats = 'on';
    case false
        stats = 'off';
end

%%

[t, X] = solver(f, [0 T], x0, odeset('RelTol',RelTol,'Stats',stats,'MaxStep',MaxStep)); %,'OutputFcn',@odeplot
n = size(X,2)/2;
x = X(:, 1:n);
v = X(:, n+1:2*n);
va = zeros(size(X));
for it = 1:length(t)
    va(it, :) = f(t(it), X(it, :).').';
end
a = va(:, n+1:2*n);

if visible
    switch output
        case 'x'
            multiplePlot(t, x, "x");
        case 'v'
            multiplePlot(t, v, "v");
        case 'a'
            multiplePlot(t, a, "a");
    end
end

switch output
    case 'x'
        xout = x;
    case 'v'
        xout = v;
    case 'a'
        xout = a;
end


    function multiplePlot(t, x, xname)
        n = size(x, 2);
        fig = figure;
        points = {};
        axes = [];
        for k=1:n
            axes(end+1) = subplot(n, 1, k);
            set(gca, 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
            hold on
            plot(t, x(:,k), 'Color', [0 0.4470 0.7410]);
            points{end+1} = plot(t, x(:,k), 'o', 'Color', [0 0.4470 0.7410], 'Visible', 'off');
            hold off
%             xlabel('t');
            ylabel(xname + k);
        end
        xlabel('t');
        linkaxes(axes,'x')
        set(fig, 'WindowButtonDownFcn', @(~, ~) tooglePoints(points));
    end

    function tooglePoints(points)
        toggle = isequal(get(points{1}, 'Visible'), 'off');
        for p = [points{:}]
            set(p, 'Visible', toggle);
        end
    end

end