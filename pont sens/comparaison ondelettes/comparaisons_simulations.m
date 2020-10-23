ti = -10;
tf = 20;
ti = -4;
tf = 10;
t1 = 0;
t2 = 3;
dt = 0.01;

masse_ajoutee = 0.2;

fex = 4; % frequence d'excitation
fp = 6; % frequence propre
zeta = 0.02; % coefficient d'amortissement
fpc = fp * (sqrt(1-zeta^2) + 1i*zeta);

Aex = [-0.5, 1 * exp(2i*pi*0.1), 0.2 * exp(2i*pi*0.7), 0.1 * exp(2i*pi*0.75)];
Ap1 = 0.8 * exp(2i*pi*0.9);
Ap2 = 1.2 * exp(2i*pi*0.5);

%%

t = ti:dt:tf;
x = Aex * exp(2i*pi*fex*[0;1;2;3]*(t-t1)) .* (t>=t1 & t<=t2);
x = x + Ap1 * exp(2i*pi*fpc/sqrt(1+masse_ajoutee^2)*(t-t1)) .* (t>=t1 & t<=t2);
x = x + Ap2 * exp(2i*pi*fpc*(t-t2)) .* (t>=t2);
x = real(x);

% x = sin(2*pi*4*t + 2*pi*rand());
% x = x .* (t<3);

figure;
plt = plot(t, x);


%%

fmin = 3;
fmax = 13;
Q = 10;
MaxRidges = 1;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MaxRidges', MaxRidges);