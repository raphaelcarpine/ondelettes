pendule = TMDpendule(9.81^2, 1, 9.81, @(theta,omega) 1*sign(omega));
pendule.reponseLibre(0, 2, 400);

%%
mr = TMDmasseressort(1, 1, @(x, v) 0.01*sign(v));
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
m0 = 100;
m1 = 1;
k0 = 100;
k1 = (1/1.01)^2;
c0 = 0;
c1 = 2/1.01*sqrt(3*0.01/8/1.01);
mr = TMDmasseressort(m1, k1, @(x, v) c1*v);
tour = Structure(m0, k0, @(x,v) 0*v, {{mr, 1}});
[t, X] = tour.reponseLibre(0, 1, 500, true);
% x = X(:, 1);
% CWT(t, x, 1/(2*pi)*exp(linspace(-0.1, 0.1, 200)), 2, 'fourier', 1);
tour.diagrammeBode(1, 1, 1/(2*pi)*logspace(-0.3, 0.3, 100), 100, true);








