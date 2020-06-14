close all
clear all

% param
% freq = @(x) 6 - abs(x);
% freq = @(x) min(5 - 0.1*log(abs(x)), 7);
freq = 6;

m_pont = 1;
m_train = 0.1;
zeta = 0.03;

T = 20;
t1 = 3;

nb_wagons = 20;

%% excitation
Cf = 400;

Lbogies = 18.70;
lbogies = 4;
essieux = Lbogies * (0:nb_wagons);
essieux = [essieux, essieux + lbogies];

v = 78.47;
L = 17.5;
deformee = @(x) (x>=0 && x<=L) .* (1-cos(2*pi*x/L))/2;
F = @(x) 0;
for k = 1:length(essieux)
    F = @(x) F(x) + deformee(x-essieux(k));
end
F = @(t) Cf * F(v*(t - t1));

% F = @(t) Cf * sin(0.42*2*pi*t);
% F = @(t) Cf * sin(4.2*2*pi*t);
% eps = 0.05;
% F = @(t) 1/eps * (mod(4.2*t+eps/2, 1) <= eps);
% F = @(~) 0;

%% eq diff


M = @(t) m_pont + m_train*(F(t)>0);
K = @(x) m_pont * (2*pi*freq).^2 *(1 + 50000*x.^2);
C = zeta * 2 * sqrt(K(0)/m_pont);

x0 = 0;
v0 = 0;

if true
    t = 0:0.01:T;
    f = nan(size(t));
    for k = 1:length(f)
        f(k) = F(t(k));
    end
    figure;
    plt0 = plot(t, f);
    xlabel('t');
    ylabel('f');
    WaveletMenu('WaveletPlot', plt0);
end


D = @(t, X) [X(2); -1/M(t) * (K(X(1))*X(1) + C*X(2) - F(t))];

[t0, X] = ode45(D, [0, T], [x0; v0]);
t = 0:0.01:T;
X = interp1(t0, X, t);
X = transpose(X);

V = nan(size(X));
for k = 1:size(X, 2)
    V(:, k) = D(t(k), X(:, k));
end

x = X(1, :);
v = X(2, :);
a = V(2, :);


figure;
plt = plot(t, x);
xlabel('t');
ylabel('x');

%% ondelette
Q = 10;
fmin = 2;
fmax = 15;

WaveletMenu('WaveletPlot', plt, 'Q', Q, 'Fmin', fmin, 'Fmax', fmax);







