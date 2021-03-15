d = 1;
L = 5.2;
N = 30;
c = 1;

d = 3;
L = 10;
N = 20;
c = 1;


k = 2;

phi_k = @(x) sin(k*pi*x/L);

sndMb = @(t) sum( phi_k((0:N-1)*d + c*t) .* (0 <= ((0:N-1)*d + c*t) & ((0:N-1)*d + c*t) <= L) );


Dt = L/c;
t1 = -(N-1)*d/c;
t2 = Dt;

T = linspace(t1-2*Dt, t2+2*Dt, 10000);
SndMb = nan(size(T));
for it = 1:length(T)
    SndMb(it) = sndMb(T(it));
end


figure('Name', sprintf('k=%d (d=%g, L=%g, N=%d, c=%g)', [k, d, L, N, c]));
plot(T, SndMb);
xlim([T(1), T(end)]);
xticks([t1, t1+Dt, t2-Dt, t2])
xticklabels({'t_1', 't_1 + \Deltat', 't_2 - \Deltat', 't_2'})
xline(t1, '--r');
xline(t1+Dt, '--r');
xline(t2-Dt, '--r');
xline(t2, '--r');
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 13);
% ylabel('$\langle\mbox{\boldmath$\varphi$}_k, \mbox{\boldmath$\mu$}_q(t)\rangle$', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('Time');
ylabel('$\langle\mbox{\boldmath$\phi$}_k,\, F_g\rangle (t) / \hat m_k$', 'Interpreter', 'latex', 'FontSize', 13);
yticks([]);



%%


t0 = 0;
tf = 10;
t1 = 1.5;
t2 = 6.5;
Dt = 0.7;

fig = figure('Name', 'réponse theorique du pont');
ax = axes(fig);
hold(ax, 'on');

T1 = linspace(t0, t1, 10);
X1 = zeros(size(T1));
l = plot(ax, T1, X1);
c = get(l, 'Color');

T2 = linspace(t1+Dt, t2-Dt, 1000);
X2 = [7+0.5*ones(size(T2)); 5+0.5*sin(10*T2+2); 3+0.2*sin(20*T2-0.5);...
    -1+0.9*sin(32*T2+0.5).*exp(-(T2-t1-Dt)/1); -3+0.8*sin(45*T2+3).*exp(-(T2-t1-Dt)/0.5);...
    -5+0.3*sin(70*T2-1).*exp(-(T2-t1-Dt)/0.4)];
plot(ax, T2, X2, 'Color', c);
xt = (t1+t2)/2 * ones(1, 7);
yt = [-6 -4 -2 0 2 4 6];
text(xt, yt, '+', 'FontSize', 12, 'HorizontalAlignment', 'center');
xt = (t1+t2)/2 * ones(1, 2);
yt = [-7 1];
text(xt, yt, '...', 'FontSize', 12, 'HorizontalAlignment', 'center');

T3 = linspace(t2, tf, 1000);
X3 = [3+0.9*sin(32*T2-5).*exp(-(T2-t1-Dt)/1); 1+0.3*sin(45*T2+0).*exp(-(T2-t1-Dt)/0.5);...
    -1+0.9*sin(70*T2-10).*exp(-(T2-t1-Dt)/0.4)];
plot(ax, T3, X3, 'Color', c);
xt = (t2+tf)/2 * ones(1, 3);
yt = [-2 0 2];
text(xt, yt, '+', 'FontSize', 12, 'HorizontalAlignment', 'center');
xt = (t2+tf)/2;
yt = -3;
text(xt, yt, '...', 'FontSize', 12, 'HorizontalAlignment', 'center');

patch([t1, t1, t1+Dt, t1+Dt], [-8 8 8 -8], c, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
patch([t2-Dt, t2-Dt, t2, t2], [-8 8 8 -8], c, 'EdgeColor', 'none', 'FaceAlpha', 0.2);


xlim([t0, tf]);
ylim([-8, 8]);
xticks([t1, t1+Dt, t2-Dt, t2])
xticklabels({'t_1', 't_1+\Deltat', 't_2-\Deltat', 't_2'})
xline(t1, '--r');
xline(t1+Dt, '--r');
xline(t2-Dt, '--r');
xline(t2, '--r');
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 13);
% ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('Time');
ylabel('Deflection');
yticks([]);








