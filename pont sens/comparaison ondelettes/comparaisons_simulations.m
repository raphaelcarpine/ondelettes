ti = -10;
tf = 20;
t1 = 2;
t2 = 5;
dt = 0.01;

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
x = x + exp(2i*pi*fpc*(t-t1)) .* (t>=t1 & t<=t2);
x = x + exp(2i*pi*fpc*(t-t2)) .* (t>=t2);
x = real(x);

figure;
plt = plot(t, x);


%%

fmin = 1;
fmax = 8;
Q = 10;
MaxRidges = 3;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MaxRidges', MaxRidges);