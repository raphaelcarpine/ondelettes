T = 200;
t = linspace(-T, T, 10000);
A = 1-0.5*sin(0.1*t);
% A = 1-0.1*t;
% A = exp(-0.2*t);
% A = max(A, 0);
% A = abs(A);
f = 1/(2*pi);
x = A .* cos(2*pi*f*t);


% plot
fig =  figure;
ax = axes(fig);
plt = plot(ax, t, x);


% ondelette%
Q = 30;
MaxRidges = 1;
MaxParallelRidges = 1;
fmin = 0.5*f;
fmax = 1.5*f;
NbFreq = 500;

ct = 3;
cf = 5;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
    , 'CtEdgeEffects', ct);


%% figure

omega = 1;

% lineaire
A = 1-0.1*t;
Q = 30;

% sinusoidale
A = 1-0.5*sin(0.1*t);

% exponentielle
A = exp(-0.2*t);
Q = 10;
"echelle log";