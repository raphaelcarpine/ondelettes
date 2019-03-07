mu = 0.01;
omega0 = 2*pi;
omega1 = 2*pi/(1+mu);
zeta0 = 0.;
Zeta1 = linspace(0, 2, 10000);

racines = zeros(4, length(Zeta1));
for k = 1:length(Zeta1)
    zeta1 = Zeta1(k);
    racines(:, k) = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
end

%%
f = figure;
ax = axes('Parent', f, 'position', [0.13 0.39  0.77 0.54]);
hold(ax, 'on');
plot(real(racines)', imag(racines)', '.', 'Parent', ax);

zeta1 = Zeta1(1);
b = uicontrol('Parent',f,'Style','slider', 'value',zeta1,...
    'min', Zeta1(1), 'max', Zeta1(end),'Position',[81,54,419,23]);
bgcolor = f.Color;
b2 = uicontrol('Parent',f,'Style','text', 'String',num2str(zeta1),...
    'BackgroundColor',bgcolor,'Position',[240,25,100,23]);

racines = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
points = plot(real(racines), imag(racines), 'o');

hold(ax, 'off');

%%
freqs = logspace(-1, 1, 10000);

bode = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs);

f2 = figure;
ax2 = axes('Parent', f2);
fbode = loglog(freqs, abs(bode)', 'Parent', ax2);

%%
b.Callback = @(es,ed) updateZeta(b2, points, mu, omega0, omega1, zeta0, es.Value, fbode, freqs);

function updateZeta(b2, points, mu, omega0, omega1, zeta0, zeta1, fbode, freqs)
racines = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
set(points, 'XData', real(racines), 'YData', imag(racines));
set(b2, 'String', num2str(zeta1));
bode = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs);
set(fbode, 'YData', abs(bode)');
end
