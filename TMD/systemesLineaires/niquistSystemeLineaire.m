mu = 1;
omega0 = 2*pi;
omega1 = 2*pi/(1+mu);
omega1 = 2*pi/2.01;
zeta0 = 0;
Zeta1 = linspace(0, 2, 100000);
Zeta1 = linspace(0, 2, 100000);
sqrt(3*mu/8/(1+mu))

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
ax = axes('Parent', f);
hold(ax, 'on');
plot(real(racines)', imag(racines)', '.', 'Parent', ax);

racines = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
points = plot(real(racines), imag(racines), 'o', 'Parent', ax);
hold(ax, 'off');

%%
freqs = logspace(-1, 1, 10000);

bode = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs);

f2 = figure;
pos = get(f, 'Position');
pos(1) = pos(1) + pos(3);
set(f2, 'Position', pos);
ax2 = axes('Parent', f2);
fbode = loglog(freqs, abs(bode)', 'Parent', ax2);
xlabel(ax2, 'f');
ylabel(ax2, '|H|');
grid(ax2, 'on');

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
bode = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs);
set(fbode, 'YData', abs(bode)');
end
