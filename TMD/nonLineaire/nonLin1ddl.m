omega0 = '2*pi';
epsilon = 0.1;
n = 1;
delta = 0.05;

d2x = '-omega0^2*x-(abs(x)<=delta)*epsilon*sign(dx)*abs(dx)^n';

T = 100;
nT = 200;

x0 = 0;
v0 = 1;


waveletplots = systemeQuelconque({'x'}, {d2x}, {'omega0', 'epsilon', 'n', 'delta'}, {omega0, epsilon, n, delta},...
    x0, v0, true, 'T', T, 'nT', nT);


%ondelette
Q = 1;
MaxParallelRidges = 1;
fmin = 0.9;
fmax = 1.1;
NbFreq = 100;


WaveletMenu('WaveletPlot', waveletplots, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);


%regression
RegressionMenu;