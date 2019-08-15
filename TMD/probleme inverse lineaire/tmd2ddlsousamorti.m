mu = 0.01;
omega0 = 2*pi;
omega1 = 2*pi/(1+mu);
zeta1 = 1 *sqrt(3*mu/8/(1+mu));



x0 = [0;0];
v0 = [1;0];




t = linspace(0, 50, 1000);



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
RegressionMenu;


%% poles (a remplir)

f1 = 1.0347;
f2 = 0.95551;
lambda1 = 0.19624;%0.18829;   
lambda2 = 0.18207;%0.17672;    

p1 = -lambda1 + 1i*2*pi*f1;
p2 = -lambda1 - 1i*2*pi*f1;
p3 = -lambda2 + 1i*2*pi*f2;
p4 = -lambda2 - 1i*2*pi*f2;



%% paramètres du problème à partir des poles

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























