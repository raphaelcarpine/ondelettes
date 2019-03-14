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
% pendule = TMDpendule(9.81^2, 1, 9.81, @(theta,omega) 1*sign(omega)*abs(omega)^(0));
% [t, x] = pendule.reponseLibre(pi-0.01, 0, 350);
pendule = TMDpendule(9.81^2, 1, 9.81, @(theta,omega) 1*omega);
[t, x] = pendule.reponseLibre(pi-0.01, 0, 1000);

Q = 2;
fmin = 0.01;
fmax = 0.2;
NbFreq = 600;

T = linspace(t(1), t(end), 1000*20/(2*pi)*2);
X = interp1(t, x, T);

WaveletMenu(fmin,fmax,NbFreq, T, X);

%%
mr = TMDmasseressort(1, 1, @(x,v) 0.1*v);
[t, x] = mr.reponseLibre(0, 1, 1000);

Q = 2;
fmin = 0.01;
fmax = 0.2;
NbFreq = 100;

T = linspace(t(1), t(end), 1000*30/(2*pi)*2);
X = interp1(t, x, T);


WaveletMenu(fmin,fmax,NbFreq, T, X);

%%
t = linspace(-10, 10, 100000);
N = -0.4;
figure;
hold on
for n=N
    plot(t, real((1i./(t+1i)).^(n+1)));
    plot(t, imag((1i./(t+1i)).^(n+1)));
end
hold off

%%
fmin = 0.5;
fmax = 1.5;
NbFreq = 50;

T = linspace(0, 100, 100*40+1);
T = T(1:(end-1));
X = (sin(2*pi*T)+0).*(T<50);

figure;
ax = plot(T, X);

WaveletMenu(fmin,fmax,NbFreq, 'WaveletPlot', ax);

%%
fmin = 0.5;
fmax = 1.5;
NbFreq = 50;

T = linspace(0, 100, 100*40+1);
T = T(1:(end-1));
X = (sin(2*pi*T)+0).*exp(-T/10);

figure;
ax = plot(T, X);

WaveletMenu(fmin,fmax,NbFreq, 'WaveletPlot', ax);


%%
mu = 0.01;
m0 = 1/mu;
m1 = 1;
k0 = 1/mu;
k1 = (1/(1+mu))^2;
c0 = 0;
c1 = 2/(1+mu)*sqrt(3*mu/8/(1+mu));
% c1 = 2/(1+mu)*sqrt(mu/(1+mu));
mr = TMDmasseressort(m1, k1, @(x, v) 0.1*sign(v)*abs(v)^2);
tour = Structure(m0, k0, @(x,v) 0*v, {{mr, 1}});
[t, X] = tour.reponseLibre(0, 1, 500, true);
x = X(:,1)';

T = linspace(t(1), t(end), 500/2/pi*200);
X = interp1(t, x, T);

fmin = 0.5/2/pi;
fmax = 1.5/2/pi;
NbFreq = 200;

WaveletMenu(fmin,fmax,NbFreq, T, X);









