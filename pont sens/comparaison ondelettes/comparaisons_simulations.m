ti = -10;
tf = 20;
ti = -2;
tf = 10;
t1 = 2;
t2 = 5;
dt = 0.01;

masse_ajoutee = 0.2;

fex = 4; % frequence d'excitation
fp = 6; % frequence propre
zeta = 0.02; % coefficient d'amortissement
fpc = fp * (sqrt(1-zeta^2) + 1i*zeta);

Aex = 1 * exp(2i*pi*0.1);
Ap1 = 0.8 * exp(2i*pi*0.9);
Ap2 = 0.8 * exp(2i*pi*0.5);

%%

t = ti:dt:tf;
x = exp(2i*pi*fex*(t-t1)) .* (t>=t1 & t<=t2);
x = x + exp(2i*pi*fpc/sqrt(1+masse_ajoutee^2)*(t-t1)) .* (t>=t1 & t<=t2);
x = x + exp(2i*pi*fpc*(t-t2)) .* (t>=t2);
x = real(x);

x = sin(2*pi*4*t + 2*pi*rand());
% x = x .* (t<3);

figure;
plt = plot(t, x);


%%

fmin = 3;
fmax = 8;
Q = 10;
MaxRidges = 1;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MaxRidges', MaxRidges);