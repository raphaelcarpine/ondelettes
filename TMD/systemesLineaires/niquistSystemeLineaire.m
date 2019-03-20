function niquistSystemeLineaire()

clear all;
close all;

plotValue = 'X';


mu = 0.01;
omega0 = 2*pi;
zeta0 = 0.0;
omega1 = 2*pi/(1+mu);
nZeta1 = 10000;
Zeta10 = 0;
Zeta11 = 0.3;
Zeta1 = linspace(Zeta10, Zeta11, nZeta1);

%diagrame de bode
freqs = logspace(-0.5, 0.5, 10000);

%reponse temporelle
T = 200;
nT = 10000;
t = linspace(0, T, nT);
x0 = [0 0];
v0 = [0 1];

%ondelette
Q = 1;
fmin = 0.9;
fmax = 1.1;
NbFreq = 200;

%%
racines = nan(4, length(Zeta1));
for k = 1:length(Zeta1)
    zeta1 = Zeta1(k);
    racines(:,k) = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
end
racines = racines(~isnan(racines));
racines = racines(imag(racines) >= 0); %on garde seulement les poles à partie imaginaire positive pour l'affichage


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
poles = plot(real(racines)', imag(racines)', '.', 'Parent', ax, 'Color', colors(index, :));
xlabel('Re');
ylabel('Im');

racines = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
racines = racines(imag(racines) >= 0);
points = plot(real(racines), imag(racines), 'o', 'Parent', ax);
hold(ax, 'off');

%reponse temporelle
ax3 = axes('Parent', f, 'Position', [0.05 0.25 0.45 0.3]);

[X, V, A, Xtmd, Vtmd, Atmd] = reponseTemporelleSystemeLineaire(x0, v0, t, mu, omega0, omega1, zeta0, zeta1);

reponseTemp = plot(t, eval(plotValue), 'Parent', ax3);
xlabel('t');
% ylabel(plotValue);


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

%plan poincaré
ax2 = axes('Parent', f, 'Position', [0.58 0.1 0.4 0.25]);
poinca = plot(X, V, 'Parent', ax2);
xlabel(ax2, 'x');
ylabel(ax2, 'v');


%%
b = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.05 0.14 0.45 0.03],'Style','slider',...
    'value',zeta1, 'min', Zeta10, 'max', Zeta11);
bgcolor = f.Color;
b2 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.31 0.09 0.15 0.05],'Style','edit',...
    'String',num2str(zeta1), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.2 0.09 0.1 0.05],'Style','text',...
    'String','zeta1', 'BackgroundColor',bgcolor);

b.Callback = @(es,ed) updateZeta(es.Value);
b2.Callback = @(es,ed) updateZeta(es.String);

bT = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.41 0.18 0.09 0.04],'Style','edit',...
    'String',num2str(T), 'BackgroundColor',bgcolor);

bplotvalue = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.01 0.18 0.09 0.04],'Style','edit',...
    'String',plotValue, 'BackgroundColor',bgcolor);

bT.Callback = @(~,~) updateT();
bplotvalue.Callback = @(~,~) updateT();

bmu = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.01 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(mu), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.01 0.05 0.09 0.04],'Style','text',...
    'String','mu', 'BackgroundColor',bgcolor);

bw = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.11 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(omega1/omega0), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.11 0.05 0.09 0.04],'Style','text',...
    'String','w1/w0', 'BackgroundColor',bgcolor);

bzeta0 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.21 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(zeta0), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.21 0.05 0.09 0.04],'Style','text',...
    'String','zeta0', 'BackgroundColor',bgcolor);

bZeta10 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.31 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(Zeta10), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.31 0.05 0.09 0.04],'Style','text',...
    'String','Zeta1(1)', 'BackgroundColor',bgcolor);

bZeta11 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.41 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(Zeta11), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.41 0.05 0.09 0.04],'Style','text',...
    'String','Zeta1(end)', 'BackgroundColor',bgcolor);

bmu.Callback = @(~,~) updateAll();
bw.Callback = @(~,~) updateAll();
bzeta0.Callback = @(~,~) updateAll();
bZeta10.Callback = @(~,~) updateAll();
bZeta11.Callback = @(~,~) updateAll();

%%
WaveletMenu('fmin',fmin,'fmax',fmax,'NbFreq',NbFreq, 'WaveletPlot', reponseTemp, 'Q', Q);

%%
    function updateT()
        T = eval(get(bT, 'String'));
        plotValue = get(bplotvalue, 'String');
        [X, V, A, Xtmd, Vtmd, Atmd] = reponseTemporelleSystemeLineaire(x0, v0, linspace(0, T, nT), mu, omega0, omega1, zeta0, zeta1);
        set(reponseTemp, 'XData', linspace(0, T, nT), 'YData', eval(plotValue));
        set(poinca, 'XData', X, 'YData', V);
    end

    function updateZeta(zeta)
        if nargin == 1
            if isnumeric(zeta)
                zeta1 = zeta;
                set(b2, 'String', num2str(zeta1));
            else
                zeta1 = eval(zeta);
                set(b, 'Value', zeta1);
            end
        else
            zeta1 = eval(get(b2, 'String'));
            set(b, 'Value', zeta1);
        end
        racines = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta1);
        racines = racines(imag(racines) >= 0);
        set(points, 'XData', real(racines), 'YData', imag(racines));
        [bode0, bode1] = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs);
        set(fbode(1), 'YData', abs(bode0)');
        set(fbode(2), 'YData', abs(bode1)');
        % set(fbode(3), 'YData', angle(bode1./bode0)');
        
        updateT();
    end

    function updateAll()
        mu = eval(get(bmu, 'String'));
        omega1 = omega0*eval(get(bw, 'String'));
        zeta0 = eval(get(bzeta0, 'String'));
        Zeta10 = eval(get(bZeta10, 'String'));
        Zeta11 = eval(get(bZeta11, 'String'));
        Zeta1 = linspace(Zeta10, Zeta11, nZeta1);
        
        set(b, 'min', Zeta10, 'max', Zeta11);
        
        racines = nan(4, nZeta1);
        for k = 1:nZeta1
            zeta = Zeta1(k);
            racines(:,k) = polesSystemeLineaire(mu, omega0, omega1, zeta0, zeta);
        end
        racines = racines(~isnan(racines));
        racines = racines(imag(racines) >= 0); %on garde seulement les poles à partie imaginaire positive pour l'affichage
        set(poles, 'XData', real(racines)', 'YData', imag(racines)');
        
        updateZeta();
    end



end
