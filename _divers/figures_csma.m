%% ondelette

Q = 5;

T = 30;
t = linspace(-T, T, 1000);
psi = 1/(Q*sqrt(pi))*exp(1i*t-t.^2/(4*Q^2));
omega = linspace(0, 2, 1000);
Tpsi = 2*exp(-Q^2*(omega-1).^2);


fig = figure;
plot(t, real(psi));
hold on
plot(t, imag(psi));
xlabel('t');
ylabel('\psi');
legend({'Re\psi', 'Im\psi'});
fig.Position(3:4) = [380 280];


fig = figure;
plot(omega, real(Tpsi));
hold on
plot(omega, imag(Tpsi));
% plot(omega, abs(Tpsi));
xlabel('\omega');
ylabel('TF[\psi]');
legend({'Re TF[\psi]', 'Im TF[\psi]'});
ylim([0 2.2]);
fig.Position(3:4) = [380 280];


%% discrétisation excitation

T = 1;
dt = 0.05;
t = 0:dt:T;

w = randn(size(t));


t = [1;1;1]*t;
t = t(:);
w = [0;1;0]*w;
w = w(:);

fig = figure;
plot(t, w);
xlabel('t');
ylabel('w');
xticks([]);
yticks([]);
xlim([dt/2, T-dt/2]);
ylim([-1 1]*max(abs(get(gca, 'YLim'))));
fig.Position(3:4) = [380 280];