clear all;
close all;


mu = 0.01;
omega0 = 2*pi;
zeta0 = 0;
omega1 = 2*pi/(1+mu);
% omega1 = 2*pi/(1+mu)*0.99;
c1 = 1;
minC1 = 0;
maxC1 = 2;
n1 = 1;
minN1 = 0;
maxN1 = 5;
f1 = @(c, n) @(x,v) c*sign(v)*abs(v)^n;

%reponse temporelle
Tf = 100;
ptsParPeriode = 200;
T = linspace(0, Tf, Tf*omega0/2/pi*ptsParPeriode);
x0 = 0;
v0 = 1;

%ondelette
Q = 50;
MaxParallelRidges = 2;
fmin = 0.9;
fmax = 1.1;
NbFreq = 200;

%%
m0 = 1/mu;
m1 = 1;
k0 = omega0^2/mu;
k1 = omega1^2;
c0 = 2*zeta0*sqrt(k0*m0);


%%
f = figure;


%reponse temporelle
ax3 = axes('Parent', f, 'Position', [0.1 0.27 0.85 0.68]);

mr = TMDmasseressort(m1, k1, f1(c1, n1));
tour = Structure(m0, k0, @(x,v) c0*v, {{mr, 1}});
[t, X] = tour.reponseLibre(x0, v0, Tf, true, 'visible', false);
x = X(:,1)';
X = interp1(t, x, T);

reponseTemp = plot(T, X, 'Parent', ax3);
xlabel('t');
ylabel('x');


%%
b = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.05 0.12 0.75 0.03],'Style','slider',...
    'value',c1, 'min', minC1, 'max', maxC1);
bgcolor = f.Color;
b3 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.8 0.12 0.05 0.05],'Style','text',...
    'String','c1 = ', 'BackgroundColor',bgcolor);
b2 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.86 0.12 0.12 0.05],'Style','edit',...
    'String',num2str(c1), 'BackgroundColor',bgcolor);

bn = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.05 0.05 0.75 0.03],'Style','slider',...
    'value',n1, 'min', minN1, 'max', maxN1);
bgcolor = f.Color;
bn3 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.8 0.05 0.05 0.05],'Style','text',...
    'String','n1 = ', 'BackgroundColor',bgcolor);
bn2 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.86 0.05 0.12 0.05],'Style','edit',...
    'String',num2str(n1), 'BackgroundColor',bgcolor);



b.Callback = @(es,ed) updateC(b, b2, bn, bn2, m0, m1, k0, k1, c0, f1, es.Value, nan, x0, v0, Tf, T, reponseTemp);
b2.Callback = @(es,ed) updateC(b, b2, bn, bn2, m0, m1, k0, k1, c0, f1, str2double(es.String), nan, x0, v0, Tf, T, reponseTemp);

bn.Callback = @(es,ed) updateC(b, b2, bn, bn2, m0, m1, k0, k1, c0, f1, nan, es.Value, x0, v0, Tf, T, reponseTemp);
bn2.Callback = @(es,ed) updateC(b, b2, bn, bn2, m0, m1, k0, k1, c0, f1, nan, str2double(es.String), x0, v0, Tf, T, reponseTemp);

%%
WaveletMenu('WaveletPlot', reponseTemp, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);

%%
function updateC(b, b2, bn, bn2, m0, m1, k0, k1, c0, f1, c1, n1, x0, v0, Tf, T, reponseTemp)
if isnan(c1)
    c1 = str2double(get(b2, 'String'));
end
set(b2, 'String', num2str(c1));
set(b, 'Value', c1);
if isnan(n1)
    n1 = str2double(get(bn2, 'String'));
end
set(bn2, 'String', num2str(n1));
set(bn, 'Value', n1);

reponseTemp.Color(4) = 0.5;
drawnow;

mr = TMDmasseressort(m1, k1, f1(c1, n1));
tour = Structure(m0, k0, @(x,v) c0*v, {{mr, 1}});
[t, X] = tour.reponseLibre(x0, v0, Tf, true, 'visible', false);
x = X(:,1)';
X = interp1(t, x, T);

set(reponseTemp, 'YData', X);
reponseTemp.Color(4) = 1;
end







