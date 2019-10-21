load('maquetteTMD/donnees/mData.mat');

mesure = mesure12;

mesure = transpose(mesure);
t = mesure(1,:);
a = mesure(2,:);

%a = bandpass(a, [3, 15], 1/mean(diff(t)));

%% t0

fig = figure;
ax = axes(fig);
plot(ax, t, a);



t0 = input('t0 = ');

delete(fig);

a = a(t>=t0);
t = t(t>=t0);

%% wavelet

fig = figure;
ax = axes(fig);
plt = plot(ax, t, a);


fmin = 6;
fmax = 7.5;
Nf = 300;
Q = 50;
maxR = 2;
maxParallel = 2;


WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', Nf, 'Q', Q, 'MaxRidges', maxR, 'MaxParallelRidges', maxParallel);


%% regressions

f1 = RegressionMenu('FigureName', 'f1', 'Equation', 'f1', 'Param', 'f1', 'Param0', 6);
f2 = RegressionMenu('FigureName', 'f2', 'Equation', 'f2', 'Param', 'f2', 'Param0', 6);

A1lambda1 = RegressionMenu('FigureName', 'lambda1', 'Equation',...
    ['A1*exp(-lambda1*(x-', num2str(t0, 10), '))'],...
    'Param', 'A1 lambda1', 'Param0', [1, 1], 'Fit', 'log(y)');
A2lambda2 = RegressionMenu('FigureName', 'lambda2', 'Equation',...
    ['A2*exp(-lambda2*(x-', num2str(t0, 10), '))'],...
    'Param', 'A2 lambda2', 'Param0', A1lambda1, 'Fit', 'log(y)');
lambda1 = A1lambda1(2);
lambda2= A2lambda2(2);

%% parametres systeme


p1 = -lambda1 + 1i*2*pi*f1;
p2 = -lambda1 - 1i*2*pi*f1;
p3 = -lambda2 + 1i*2*pi*f2;
p4 = -lambda2 - 1i*2*pi*f2;


a0 = real (p1*p2*p3*p4);
a1 = -real ( p2*p3*p4 + p1*p3*p4 + p1*p2*p4 + p1*p2*p3);
a2 = real (p1*p2 + p1*p3 + p1*p4 + p2*p3 + p2*p4 + p3*p4);
a3 = -real (p1 + p2 + p3 + p4);


w1 = sqrt(a2-a0*a3/a1);
w2 = sqrt(a0)/w1;
f1 = w1/2/pi;
f2 = w2/2/pi;
zeta2 = a1 / w1^2 / (2*w2);
mu = a3/a1*w1^2 - 1;

disp(['f1 = ', num2str(f1)]);
disp(['f2 = ', num2str(f2)]);
disp(['zeta2 = ', num2str(zeta2)]);
disp(['mu = ', num2str(mu)]);
















