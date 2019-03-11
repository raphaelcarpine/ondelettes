pendule = TMDpendule(9.81^2, 1, 9.81, @(theta,omega) 1*sign(omega)*abs(omega)^(0));
pendule.reponseLibre(pi-0.05, 0, 400);

%%
mr = TMDmasseressort(1, 0, @(x, v) 0.01*sign(v)*abs(v)^(1) + sign(x));
mr.reponseLibre(0, 2, 500);

%%
pendule = TMDpendule(9.81^2, 1, 9.81, @(theta,omega) 50*omega*(abs(theta)<0.1));
tour = Structure(50, 50, @(x,v) 0.0*v, {{pendule, 1}});
tour.reponseLibreAvecTMD(0, 1, 200);

%%
l = 1;
pendule = TMDpendule(l^2, 1, l, @(theta,omega) 0);
tour = Structure(50, 50, @(x,v) 10*v, {{pendule, 1}});
tour.reponseForceeAvecTMD(0, 0, 100, @(t) [1 0], true);

%%
k = 1/1.01^2;
c = 2*sqrt(k)*sqrt(3*0.01/(8*(1+0.01)));
% c = 0.05;
% c = 10;
mr = TMDmasseressort(1, k, @(x, v) c*v);
% mr.reponseLibre(0, 1, 200);
tour = Structure(100, 100, @(x,v) 10*v, {{mr, 1}});
tour = Structure(100, 100, @(x,v) 0*v, {{mr, 1}});
[t, X] = tour.reponseLibre(1, 0, 500, true);
x = X(:, 1);
CWT(t, x, 1/(2*pi)*exp(linspace(-0.1, 0.1, 200)), 2, 'fourier', 1);
% tour.diagrammeBode(1, 1, 1/(2*pi)*exp(linspace(-0.3, 0.3, 100)), 500, true);

%%
mr = TMDmasseressort(1, 2, @(x, v) 0.01*sign(v));
[t, x] = mr.reponseLibre(0, 1, 200);
CWT(t, x, 1/(2*pi)*exp(linspace(-1, 1, 100)), 500);

%%
t = linspace(0, 100, 5000);
x = sin(2*pi*t) + 2*sin(4*pi*t+1);
CWT(t, x, logspace(-0.5, 0.5, 100), 200);

%%
mu = 0.05;
m0 = 1/mu;
m1 = 1;
k0 = 1/mu;
k1 = (1/(1+mu))^2;
c0 = 0;
c1 = 2/(1+mu)*sqrt(3*mu/8/(1+mu));
c1 = 2/(1+mu)*sqrt(mu/(1+mu));
c1 = 2/(1+mu)*sqrt(3*mu/8/(1+mu))*60227/43301;
mr = TMDmasseressort(m1, k1, @(x, v) c1*v);
tour = Structure(m0, k0, @(x,v) 0*v, {{mr, 1}});
[t, X] = tour.reponseLibre(0, 1, 100, true);
% x = X(:, 1);
% CWT(t, x, 1/(2*pi)*exp(linspace(-0.1, 0.1, 200)), 2, 'fourier', 1);
tour.diagrammeBode(1, 1, 1/(2*pi)*logspace(-0.3, 0.3, 100), 200, true);


%%
mr = TMDmasseressort(1, 1, @(x, v) 0.2*sign(v)*abs(v)^(1));
tour = Structure(50, 50, @(x,v) 0.0*v, {{mr, 1}});
tour.reponseLibre(1, 0, 200, true);

mr = TMDmasseressort(1, 0, @(x, v) 2*sign(v)*abs(v)^(1));
tour = Structure(50, 50, @(x,v) 0.0*v, {{mr, 1}});
tour.reponseLibre(1, 0, 200, true);

%%
pendule = TMDpendule(9.81^2, 1, 9.81, @(theta,omega) 1*sign(omega)*abs(omega)^(0));
[t, x] = pendule.reponseLibre(pi-0.05, 0, 350);

Q = 2;
fmin = 0.01;
fmax = 0.2;
NbFreq = 200;
WvltFreq = linspace(fmin, fmax, NbFreq);
T = linspace(t(1), t(end), 350*100/(2*pi)*2);
X = interp1(t, x, T);
WvltPlot(T, X, WvltFreq, Q);

ridge = RidgeExtract(T, X, Q, fmin, fmax, NbFreq);
figure;
plot(ridge.time{1}, abs(ridge.val{1}));
figure;
plot(abs(ridge.val{1}), ridge.freq{1});










