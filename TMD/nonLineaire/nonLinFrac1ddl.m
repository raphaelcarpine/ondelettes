omega0 = '2*pi';
epsilon = 0.1;
alpha = 1;

d2x = '-omega0^2*x-epsilon*d(alpha)x';

T = 100;
nT = 200;

x0 = 0;
v0 = 1;


waveletplots = systemeQuelconque({'x'}, {d2x}, {'omega0', 'epsilon', 'alpha'}, {omega0, epsilon, alpha},...
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