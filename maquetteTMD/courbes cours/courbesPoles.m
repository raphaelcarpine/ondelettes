%% premier pole

pa = -0.5 + 20i;
phi = 1;



t = linspace(0, 10, 1000);
X = exp(pa*t+1i*phi);
x = real(X);
A = abs(X);


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(ax, t, x);
plot(ax, t, A);
plot(ax, [t(1), t(end)], [0, 0], 'black');
plot(ax, [0, -1/real(pa)], [1, 0], 'r--', 'LineWidth', 2);
plot(ax, -1/real(pa), 0, 'r+', 'LineWidth', 2);
text(ax, -1/real(pa), -0.1, '-1/Re(p_a)', 'FontSize', 15, 'Color', 'red');
ylim(ax, [-0.5, 1]);
xlabel(ax, 't', 'FontSize', 15);
ylabel(ax, 'e^{p_at}', 'FontSize', 15);




%% premier pole

pb = -0.2 + 23i;
phi = 1;



X = exp(pb*t+1i*phi);
x = real(X);
A = abs(X);


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(ax, t, x);
plot(ax, t, A);
plot(ax, [t(1), t(end)], [0, 0], 'black');
plot(ax, [0, -1/real(pb)], [1, 0], 'r--', 'LineWidth', 2);
plot(ax, -1/real(pb), 0, 'r+', 'LineWidth', 2);
text(ax, -1/real(pb), -0.1, '-1/Re(p_b)', 'FontSize', 15, 'Color', 'red');
ylim(ax, [-0.5, 1]);
xlabel(ax, 't', 'FontSize', 15);
ylabel(ax, 'e^{p_bt}', 'FontSize', 15);


%% somme

S = exp(pa*t+1i*phi) + 1.3*exp(pb*t+1i*phi);


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(ax, t, real(S));
xlabel(ax, 't', 'FontSize', 15);
ylabel(ax, 'C_ae^{p_at} + C_be^{p_bt}', 'FontSize', 15);

