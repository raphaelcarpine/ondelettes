%% figures autocorrélation bruit blanc

N = 10000;
Dt = 0.05;
sigma = 1;

%% bruit
n = 0:N-1;
w = sigma * randn(1, N);

figure;
plot(n, w);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$\tilde w$', 'interpreter', 'latex', 'FontSize', 13);
f = gcf; f.Position(3:4) = 7/10 * [560 420];

%% autocorr
n = -N+1:N-1;
Rw = Dt * xcorr(w);

figure;
plot(n, Rw);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$R_{\tilde w}$', 'interpreter', 'latex', 'FontSize', 13);
f = gcf; f.Position(3:4) = 7/10 * [560 420];