mu = 0.001;
omega0 = 2*pi;
omega1 = 2*pi/(1+mu);
zeta0 = 0;
Zeta1 = linspace(0, 0.05, 100000);

disp(['equal pic design : zeta = ' num2str(sqrt(3*mu/8/(1+mu)))]);
disp(['amortissement max : zeta = ' num2str(sqrt(mu/(1+mu)))]);


freqs = logspace(-0.5, 0.5, 10000);

t = linspace(0, 2000, 100000);
x0 = 0;
v0 = 1;

%%
racines = nan(4, length(Zeta1));
for k = 1:length(Zeta1)
    zeta1 = Zeta1(k);
    r = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
    for i = 1:length(r)
        racines(i, k) = r(i);
    end
end

%%
f = figure;
zeta1 = Zeta1(1);
b = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.05 0.1 0.45 0.05],'Style','slider',...
    'value',zeta1, 'min', Zeta1(1), 'max', Zeta1(end));
bgcolor = f.Color;
b2 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.05 0.05 0.45 0.05],'Style','edit',...
    'String',num2str(zeta1), 'BackgroundColor',bgcolor);

%%
ax = axes('Parent', f, 'Position', [0.05 0.65 0.45 0.3]);
hold(ax, 'on');
colors = get(ax, 'ColorOrder');
index = get(ax,'ColorOrderIndex');
plot(real(racines)', imag(racines)', '.', 'Parent', ax, 'Color', colors(index, :));
xlabel('Re');
ylabel('Im');

racines = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
points = plot(real(racines), imag(racines), 'o', 'Parent', ax);
hold(ax, 'off');

%%
ax3 = axes('Parent', f, 'Position', [0.05 0.25 0.45 0.3]);

Xt = reponseTemporelleSystemeLineaire(x0, v0, t, mu, omega0, omega1, zeta0, zeta1);

reponseTemp = plot(t, Xt, 'Parent', ax3);
xlabel('t');
ylabel('x');


%%

[bode0, bode1] = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs);

axes1 = [];
fbode = [];

axes1(1) = axes('Parent', f, 'Position', [0.58 0.7 0.4 0.25]);
fbode(1) = loglog(freqs, abs(bode0)', 'Parent', axes1(end));
set(axes1(1), 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
ylabel('|H0|');

axes1(2) = axes('Parent', f, 'Position', [0.58 0.4 0.4 0.25]);
fbode(2) = loglog(freqs, abs(bode1)', 'Parent', axes1(2));
set(axes1(2), 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
ylabel('|H1|');

axes1(3) = axes('Parent', f, 'Position', [0.58 0.1 0.4 0.25]);
fbode(3) = semilogx(freqs, angle(bode1./bode0)', 'Parent', axes1(3));
set(axes1(3), 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
ylim(axes1(3), [-pi pi]);
ylabel('arg(H1/H0)');

xlabel('f');
linkaxes(axes1,'x')


%%
b.Callback = @(es,ed) updateZeta(b, b2, points, mu, omega0, omega1, zeta0, es.Value, fbode, freqs, x0, v0, t, reponseTemp);
b2.Callback = @(es,ed) updateZeta(b, b2, points, mu, omega0, omega1, zeta0, str2double(es.String), fbode, freqs, x0, v0, t, reponseTemp);


%%
function updateZeta(b, b2, points, mu, omega0, omega1, zeta0, zeta1, fbode, freqs, x0, v0, t, reponseTemp)
racines = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
set(points, 'XData', real(racines), 'YData', imag(racines));
set(b2, 'String', num2str(zeta1));
set(b, 'Value', zeta1);
[bode0, bode1] = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs);
set(fbode(1), 'YData', abs(bode0)');
set(fbode(2), 'YData', abs(bode1)');
set(fbode(3), 'YData', angle(bode1./bode0)');
Xt = reponseTemporelleSystemeLineaire(x0, v0, t, mu, omega0, omega1, zeta0, zeta1);
set(reponseTemp, 'YData', Xt);
end
