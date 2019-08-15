mu = 0.01;
omega0 = 2*pi;
omega1 = 2*pi/(1+mu);
zeta1 = 2 *sqrt(3*mu/8/(1+mu));




x0 = [0;0];
v0 = [1;0];




t = linspace(0, 50, 1000);

poles = polesSystemeLineaire(mu, omega0, omega1, 0, zeta1);



%% reponse temporelle et ondelette

[x, v, a, xtmd, vtmd, atmd] = reponseTemporelleSystemeLineaire(x0, v0, t, mu, omega0, omega1, 0, zeta1);

fig = figure;
ax = axes(fig);
plt = plot(t, x, 'Parent', ax);



%ondelette
Q = 25;
MaxParallelRidges = 1;
fmin = 0.9;
fmax = 1.1;
NbFreq = 100;


WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);


%regression
Equation = 'sqrt (a1^2*exp(-2*l1*x) + a2^2*exp(-2*l2*x) + 2*a1*a2*exp(-(l1+l2)*x)*cos(dom*x+d))';
Param = 'a1 a2 l1 l2 dom d';
Param0 = '0.1 0.1 0.1 0.1 0.1 0';
Fit = 'log(y)';
RegressionMenu('Equation', Equation, 'Param', Param, 'Param0', Param0, 'Fit', Fit);


%% paramètres obtenus avec l'amplitude (a remplir)

a1 = 0.058923;
a2 = 0.21702;
l1 = 0.59098;
l2 = 0.16183;
dOm = 0.027597;
D = -3.0514;



%% extraction avec la fréquence 

[tf, freqr] = PlotExtract();


dpsi = @(dOm, D, t) dOm * (1 + a1/a2*exp((l2-l1)*t) .* (cos(dOm*t+D) - (l2-l1)/dOm*sin(dOm*t+D)))...
    ./ ( (a1/a2*exp((l2-l1)*t) + cos(dOm*t+D)).^2 + sin(dOm*t+D).^2);


%%
%premier cas
dpsit = dpsi(dOm, D, tf);

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tf, 2*pi*freqr, 'Parent', ax);
plot(tf, 2*pi*freqr-dpsit, 'Parent', ax);

%deuxieme cas
dpsit = dpsi(-dOm, -D, tf);

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tf, 2*pi*freqr, 'Parent', ax);
plot(tf, 2*pi*freqr-dpsit, 'Parent', ax);

%%
dOm = -dOm;
D = -D;

dpsit = dpsi(dOm, D, tf);

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tf, 2*pi*freqr, 'Parent', ax);
plot(tf, 2*pi*freqr-dpsit, 'Parent', ax);



RegressionMenu;

%%
Om1 = 6.2773;
Om2 = Om1 + dOm;


%% paramètres du problème à partir des poles

p1 = -l1 + 1i*Om1;
p2 = -l1 - 1i*Om1;
p3 = -l2 + 1i*Om2;
p4 = -l2 - 1i*Om2;


a0 = real (p1*p2*p3*p4);
a1 = -real ( p2*p3*p4 + p1*p3*p4 + p1*p2*p4 + p1*p2*p3);
a2 = real (p1*p2 + p1*p3 + p1*p4 + p2*p3 + p2*p4 + p3*p4);
a3 = -real (p1 + p2 + p3 + p4);


omega02 = sqrt(a2-a0*a3/a1)
(omega02-omega0)/omega0

omega12 = sqrt(a0)/omega02
(omega12-omega1)/omega1

zeta12 = a1 / omega02^2 / (2*omega12)
(zeta12-zeta1)/zeta1

mu2 = a3/a1*omega02^2 - 1
(mu2-mu)/mu























