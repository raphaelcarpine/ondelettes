omega0 = '2*pi';
epsilon = 0.1;
n = 1;
omegaK = '2*pi';

d2x = ' -omega0^2*x - omegaK^2*(x-z)';
dz = 'sign(x-z) * abs(omegaK^2/epsilon*(x-z))^(1/n)';

T = 100;
nT = 200;

x0 = 0;
v0 = 1;
z0 = 0;


waveletplots = systemeQuelconque({'x'}, {d2x}, {'omega0', 'epsilon', 'n', 'omegaK'}, {omega0, epsilon, n, omegaK},...
    x0, v0, false, 'variablesNonInertielles', {'z'}, 'equationsNonInertielles', {dz}, 'X0NonInertiel', z0,...
    'T', T, 'nT', nT);


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