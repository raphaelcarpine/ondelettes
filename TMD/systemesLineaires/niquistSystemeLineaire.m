clear all;
close all;


mu = 0.01;
omega0 = 2*pi;
zeta0 = 0.0;
omega1 = 2*pi/(1+mu);
Zeta1 = linspace(0, 0.3, 100000);

disp(['equal pic design : zeta = ' num2str(sqrt(3*mu/8/(1+mu)))]);
disp(['amortissement max : zeta = ' num2str(sqrt(mu/(1+mu)))]);

%diagrame de bode
freqs = logspace(-0.5, 0.5, 10000);

%reponse temporelle
t = linspace(0, 200, 10000);
x0 = 0;
v0 = 1;

%ondelette
Q = 1;
fmin = 0.9;
fmax = 1.1;
NbFreq = 200;

%%
racines = nan(4, length(Zeta1));
for k = 1:length(Zeta1)
    zeta1 = Zeta1(k);
    r = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
    r = r(imag(r) >= 0); %on garde seulement les poles � partie imaginaire positive pour l'affichage
    for i = 1:length(r)
        racines(i, k) = r(i);
    end
end

%%
f = figure;
f.Position(2) = f.Position(2) - 0.5*f.Position(4);
f.Position(3:4) = f.Position(3:4) + 0.5*f.Position(3:4);


zeta1 = Zeta1(1);

%poles
ax = axes('Parent', f, 'Position', [0.05 0.65 0.45 0.3]);
hold(ax, 'on');
colors = get(ax, 'ColorOrder');
index = get(ax,'ColorOrderIndex');
plot(real(racines)', imag(racines)', '.', 'Parent', ax, 'Color', colors(index, :));
xlabel('Re');
ylabel('Im');

racines = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
racines = racines(imag(racines) >= 0);
points = plot(real(racines), imag(racines), 'o', 'Parent', ax);
hold(ax, 'off');

%reponse temporelle
ax3 = axes('Parent', f, 'Position', [0.05 0.25 0.45 0.3]);

[Xt, Vt] = reponseTemporelleSystemeLineaire(x0, v0, t, mu, omega0, omega1, zeta0, zeta1);

reponseTemp = plot(t, Xt, 'Parent', ax3);
xlabel('t');
ylabel('x');


%bode
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

% axes1(3) = axes('Parent', f, 'Position', [0.58 0.1 0.4 0.25]);
% fbode(3) = semilogx(freqs, angle(bode1./bode0)', 'Parent', axes1(3));
% set(axes1(3), 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
% ylim(axes1(3), [-pi pi]);
% ylabel('arg(H1/H0)');

xlabel('f');
linkaxes(axes1,'x')

%plan poincar�
ax2 = axes('Parent', f, 'Position', [0.58 0.1 0.4 0.25]);
poinca = plot(Xt, Vt, 'Parent', ax2);
xlabel(ax2, 'x');
ylabel(ax2, 'v');


%%
b = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.05 0.14 0.45 0.03],'Style','slider',...
    'value',zeta1, 'min', Zeta1(1), 'max', Zeta1(end));
bgcolor = f.Color;
b2 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.31 0.09 0.15 0.05],'Style','edit',...
    'String',num2str(zeta1), 'BackgroundColor',bgcolor);
b3 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.2 0.09 0.1 0.05],'Style','text',...
    'String','zeta1', 'BackgroundColor',bgcolor);

b.Callback = @(es,ed) updateZeta(b, b2, points, mu, omega0, omega1, zeta0, es.Value, fbode, freqs, x0, v0, t, reponseTemp, poinca);
b2.Callback = @(es,ed) updateZeta(b, b2, points, mu, omega0, omega1, zeta0, str2double(es.String), fbode, freqs, x0, v0, t, reponseTemp, poinca);

%%
WaveletMenu('fmin',fmin,'fmax',fmax,'NbFreq',NbFreq, 'WaveletPlot', reponseTemp, 'Q', Q);

%%
function updateZeta(b, b2, points, mu, omega0, omega1, zeta0, zeta1, fbode, freqs, x0, v0, t, reponseTemp, poinca)
racines = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
racines = racines(imag(racines) >= 0);
set(points, 'XData', real(racines), 'YData', imag(racines));
set(b2, 'String', num2str(zeta1));
set(b, 'Value', zeta1);
[bode0, bode1] = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs);
set(fbode(1), 'YData', abs(bode0)');
set(fbode(2), 'YData', abs(bode1)');
% set(fbode(3), 'YData', angle(bode1./bode0)');
[Xt, Vt] = reponseTemporelleSystemeLineaire(x0, v0, t, mu, omega0, omega1, zeta0, zeta1);
set(reponseTemp, 'YData', Xt);
set(poinca, 'XData', Xt, 'YData', Vt);
end







