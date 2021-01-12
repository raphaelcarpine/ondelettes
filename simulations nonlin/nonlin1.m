%% system
m = 1;
k1 = (10*2*pi)^2;
k2 = k1/2;
x0 = 1;
c = 0.1;


kx = @(x) k1*x*(x<=x0) + (k2*(x-x0) + k1*x0)*(x>x0);

%% time
T = 50;
fe = 500;

dt = 1/fe;
t = 0:dt:T;

%% excitation

% dirac
f = zeros(size(t));
% f(1) = 1 / dt;

% noise

% noise + diracs


%% ode
x0 = 3;
v0 = 0;

X0 = [x0; v0];

D = @(t, X) [X(2); 1/m*(-kx(X(1)) - c*X(2) - f(floor(t/dt)+1))];

[~, X] = ode45(D, t, X0);
x = X(:, 1).';

%% plot
figure;
plt = plot(t, x);
xlabel('Time [s]');
ylabel('Displacement [m]');

%% wavelet
fmin = 5;
fmax = 50;
Q = 10;
