f0 = 1;
z = 0.02;
phi = [1; 0.6];
phi = [1; 0.6; -0.8];
% phi = [1; 0.6*exp(0.4i); -0.8*exp(-0.2i)];
% phi = [1; 0.6*exp(0.4i); -0.8*exp(0.2i); -1.1*exp(-0.1i)];

w0 = 2*pi*f0;

%%

dt = 0.001;
T = 15;%10;
t = 0:dt:T+1;
n = length(t);
f = 0.5*randn(size(t));
f(1) = f(1) + 1/dt;

H = @(w) 1./(w0^2 + 2i*z*w0*w - w.^2);
hilb = [ones(1, floor(n/2)+1), zeros(1, ceil(n/2)-1)];
W = 1/(n*dt)*[0:floor(n/2), -ceil(n/2)+1:-1];
X = real(ifft(phi * (H(W).*fft(f).*hilb), [], 2));

figure;
plt = plot(t, X);
xlabel('Temps [s]');
ylabel('Déplacement');
xlim([0 T]);
% ylim(max(X, [], 'all')*[-1.7 1.7]);
ylim(max(X, [], 'all')*[-1.3 1.3]);
% xticks([0 T]);
% xticklabels({'0', 'T'});
% xticks([]);
yticks([]);
ax = gca;
fig = gcf;
fig.Position(3:4) = [450 350];
hold on

ti = 0.22*T; tf = 0.3*T;
t0 = t(t >= ti & t <= tf);
X0 = X(:, t >= ti & t <= tf);
x0lim = 1.2*max(abs(X0), [], 'all');
plot([ti tf tf ti ti], x0lim*[1 1 -1 -1 1], 'k');

ax2Pos = [0.25 0.15 0.3 0.25];
ax2Pos = [0.5 0.15 0.3 0.25];
posx = @(x) ax.XLim(1) + (x-ax.InnerPosition(1))/ax.InnerPosition(3) * diff(ax.XLim);
posy = @(y) ax.YLim(1) + (y-ax.InnerPosition(2))/ax.InnerPosition(4) * diff(ax.YLim);

plot([ti, posx(ax2Pos(1))], [x0lim, posy(ax2Pos(2)+ax2Pos(4))], '--', 'Color', 0.4*[1 1 1]);
plot([ti, posx(ax2Pos(1))], [-x0lim, posy(ax2Pos(2))], '--', 'Color', 0.4*[1 1 1]);
plot([tf, posx(ax2Pos(1)+ax2Pos(3))], [x0lim, posy(ax2Pos(2)+ax2Pos(4))], '--', 'Color', 0.4*[1 1 1]);
plot([tf, posx(ax2Pos(1)+ax2Pos(3))], [-x0lim, posy(ax2Pos(2))], '--', 'Color', 0.4*[1 1 1]);

% legend({'u_1', 'u_2', 'u_3'});
legend({'u_1', 'u_2'});

ax2 = axes(fig, 'Position', ax2Pos, 'Box', true);
plot(t0, X0);
xlim([ti tf]);
ylim([-x0lim, x0lim]);
xticks([]);
yticks([]);

WaveletMenu('WaveletPlot', plt, 'fmin', 5, 'fmax', 8, 'XLim', [0 T]);

%%

phi = phi / sqrt(phi.' * phi);
phi = exp(2i)*phi;
phi = phi.^2;
Z = [zeros(size(phi)), phi];


figure;
for k = 1:length(phi)
    polarplot(angle(Z(k, :)), abs(Z(k, :)), 'LineWidth', 2, 'HandleVisibility', 'off');
    hold on
    ax = gca;
    ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
    polarplot(angle(Z(k, 2)), abs(Z(k, 2)), 'o', 'LineWidth', 2, 'MarkerSIze', 7,...
        'DisplayName', ['\phi^{(k)}_', num2str(k)]);
end

fig = gcf;
fig.Position(3:4) = [250 150];
fig.Position(3:4) = [300 300];

ax.ThetaTickLabel = {};
ax.RTickLabel = {};

% legend('FontSize', 9);




