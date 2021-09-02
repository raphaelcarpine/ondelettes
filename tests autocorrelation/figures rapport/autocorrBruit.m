%% figures autocorrélation bruit blanc

N = 10000;
Dt = 0.05;
sigma = 1;

%% bruit
n = 0:N-1;
w = sigma * randn(1, N);

figure;
plot(n, w);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 15);
ylabel('$\tilde w$', 'interpreter', 'latex', 'FontSize', 15);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ax = gca;
ax.YLim = max(abs(ax.YLim)) * [-1 1];

%% autocorr
n = -N+1:N-1;
Rw = Dt * xcorr(w);

figure;
plot(n, Rw);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 15);
ylabel('$R_{\tilde w}$', 'interpreter', 'latex', 'FontSize', 15);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ax = gca;
YLimR = ax.YLim;

ERw = zeros(size(Rw));
ERw((end-1)/2+1) = N*Dt*sigma^2;
figure;
plot(n, ERw);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 15);
ylabel('$N\Delta t\, \sigma^2 \delta_n$', 'interpreter', 'latex', 'FontSize', 15);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ax = gca;
ax.YLim = YLimR;

figure;
plot(n, Rw - ERw);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 15);
ylabel('$W$', 'interpreter', 'latex', 'FontSize', 15);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ax = gca;
ax.YLim = YLimR;

