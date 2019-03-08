mu = 0.01;
omega0 = 2*pi;
omega1 = 2*pi/(1+mu);
zeta0 = 0;
Zeta1 = linspace(0, 0.1, 100000);

disp(['equal pic design : zeta = ' num2str(sqrt(3*mu/8/(1+mu)))]);
disp(['amortissement max : zeta = ' num2str(sqrt(mu/(1+mu)))]);

racines = nan(4, length(Zeta1));
for k = 1:length(Zeta1)
    zeta1 = Zeta1(k);
    r = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
    for i = 1:length(r)
        racines(i, k) = r(i);
    end
end

%%
fb = figure;
zeta1 = Zeta1(1);
b = uicontrol('Parent',fb,'Style','slider', 'value',zeta1,...
    'min', Zeta1(1), 'max', Zeta1(end),'Position',[81,54,419,23]);
bgcolor = fb.Color;
b2 = uicontrol('Parent',fb,'Style','text', 'String',num2str(zeta1),...
    'BackgroundColor',bgcolor,'Position',[240,25,100,23]);

%%
f = figure;
pos = get(f, 'Position');
pos(1) = 200;
set(f, 'Position', pos);
ax = axes('Parent', f);
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
freqs = logspace(-0.5, 0.5, 10000);

[bode0, bode1] = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs);

f2 = figure;
pos = get(f, 'Position');
pos(1) = pos(1) + pos(3) + 5;
pos(2) = pos(2) - pos(4);
pos(4) = pos(4)*2;
set(f2, 'Position', pos);

axes1 = [];
fbode = [];

axes1(end+1) = subplot(3, 1, 1);
fbode(end+1) = loglog(freqs, abs(bode0)', 'Parent', axes1(end));
set(axes1(end), 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
ylabel('|H0|');

axes1(end+1) = subplot(3, 1, 2);
fbode(end+1) = loglog(freqs, abs(bode1)', 'Parent', axes1(end));
set(axes1(end), 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
ylabel('|H1|');

axes1(end+1) = subplot(3, 1, 3);
fbode(end+1) = semilogx(freqs, angle(bode1./bode0)', 'Parent', axes1(end));
set(axes1(end), 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
ylim(axes1(end), [-pi pi]);
ylabel('arg(H1/H0)');

xlabel('f');
linkaxes(axes1,'x')


%%
pos = get(f, 'Position');
pos(4) = 100;
pos(2) = pos(2) - pos(4) - 100;
set(fb, 'Position', pos); 

b.Callback = @(es,ed) updateZeta(b2, points, mu, omega0, omega1, zeta0, es.Value, fbode, freqs);

%%
set(f, 'CloseRequestFcn', closeFunction(f, f2, fb));
set(f2, 'CloseRequestFcn', closeFunction(f, f2, fb));
set(fb, 'CloseRequestFcn', closeFunction(f, f2, fb));

function func = closeFunction(f, f2, fb)
    function Func(~, ~)
        delete(f);
        delete(f2);
        delete(fb);
    end
func = @Func;
end


%%
function updateZeta(b2, points, mu, omega0, omega1, zeta0, zeta1, fbode, freqs)
racines = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
set(points, 'XData', real(racines), 'YData', imag(racines));
set(b2, 'String', num2str(zeta1));
[bode0, bode1] = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs);
set(fbode(1), 'YData', abs(bode0)');
set(fbode(2), 'YData', abs(bode1)');
set(fbode(3), 'YData', angle(bode1./bode0)');
end
