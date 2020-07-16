% clear all
% close all


%% système

m = 1;
k = @(x) (2*pi* max(1 + 0.1*log(abs(x)), 1)).^2;
c = 0.01 * 2*sqrt(k(0)/m);

phi = [1; -2; 0.5];
phi = 1;

nbDDL = 3;



%% excitation

T = 10000;
dt = 0.01;
fe = 1/dt;

t = 0:dt:T;
nt = length(t);

excitation = 'gaussien'; % 'bruit' 'dirac' 'gaussien'

if isequal(excitation, 'bruit')
    f = exp(2i*pi*rand(1, nt)); % bruit
    f(1) = 1;
    for k = 2:nt-1
        f(end+2-k) = conj(f(k));
    end
    f = ifft(f);
elseif isequal(excitation, 'dirac')
    f = 2e4*ones(1, nt); % dirac
    f = ifft(f);
elseif isequal(excitation, 'gaussien')
%     f = normrnd(0, 1, 1, nt);
%     f = - sqrt(2) *  erfcinv(2*rand(1, nt));
    f = 1e1*randn(1, nt) / sqrt(dt);
else
    error('');
end

% figure;
% histogram(f);

F = @(t) f(floor(t/dt) + 1);


%%


D = @(t, X) [X(2); -1/m * (k(X(1))*X(1) + c*X(2) - F(t))];

x0 = 0; v0 = 0;
options = odeset('MaxStep', dt);
[t0, X] = ode45(D, [0, T], [x0; v0], options);

X = interp1(t0, X, t);
X = transpose(X);

% V = nan(size(X));
% for k = 1:size(X, 2)
%     V(:, k) = D(t(k), X(:, k));
% end

x = X(1, :);
x = phi*x;

figure;
plt = plot(t, x);
xlabel('t');
ylabel('x');


%% reponse

% ondelette
Q = 10;
MaxRidges = 1;
MaxParallelRidges = 1;
fmin = 0.2;
fmax = 2;
MultiSignalMode = false;
AutocorrelationMode = true;

ct = 3;
cf = 5;

% plots

if false % f plot
    fig = figure('Name', 'excitation');
    ax = axes(fig);
    pltf = plot(ax, t, f);
    xlabel(ax, 't');
    ylabel(ax, 'x');
    
%     WaveletMenu('WaveletPlot', pltf, 'fmin', fmin, 'fmax', fmax,...
%         'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
%         , 'CtEdgeEffects', ct);
end

% x = reponseSystNddl(t, f, nbDDL, w0, zeta, C);
% ddlF = 2;
% x = syst.response(f, dt, ddlF);


WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
    , 'CtEdgeEffects', ct, 'MultiSignalMode', MultiSignalMode,...
    'AutocorrelationMode', AutocorrelationMode, 'AutocorrelationMaxLag', 200);





