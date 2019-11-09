%% simulation

f0 = 6.653;
f1 = 6.617;
zeta1 = 0.064;
mu = 0.01472;



w0 = 2*pi*f0;
w1 = 2*pi*f1;
c1 = 2*zeta1*w1;


D = @(t, XV) [XV(3); XV(4);...
    -w0^2*XV(1) - mu*w1^2*(XV(1)-XV(2)) - mu*c1*(XV(3)-XV(4));...
    -w1^2*(XV(2)-XV(1)) - c1*(XV(4)-XV(3))];

T = 10;
XV0 = [0;0;1;0];

opts = odeset('MaxStep',1e-2);

[t, XV] = ode45 (D, [0 T], XV0, opts);

t0 = t';
a = XV(:,1)';

t = linspace(0, T, length(t0));
a = interp1(t0, a, t);

% wavelet
fig = figure;
ax = axes(fig);
plt = plot(ax, t, a);


fmin = 6;
fmax = 8;
Nf = 300;
Q = 30;
maxR = 2;
maxParallel = 2;
MinModu = 10;


WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq',...
    Nf, 'Q', Q, 'MaxRidges', maxR, 'MaxParallelRidges', maxParallel, 'RidgeMinModu', MinModu);

% regression
fa = RegressionMenu('FigureName', 'fa', 'Equation', 'fa', 'Param', 'fa', 'Param0', 6);
fb = RegressionMenu('FigureName', 'fb', 'Equation', 'fb', 'Param', 'fb', 'Param0', 6);

Aala = RegressionMenu('FigureName', 'la', 'Equation', 'Aa*exp(-la*x)',...
    'Param', 'Aa la', 'Param0', [1, 1], 'Fit', 'log(y)');
Ablb = RegressionMenu('FigureName', 'lb', 'Equation', 'Ab*exp(-lb*x)',...
    'Param', 'Ab lb', 'Param0', [1, 1], 'Fit', 'log(y)');
la = Aala(2);
lb= Ablb(2);

% parametres systeme
p1 = -la + 1i*2*pi*fa;
p2 = -la - 1i*2*pi*fa;
p3 = -lb + 1i*2*pi*fb;
p4 = -lb - 1i*2*pi*fb;


a0 = real (p1*p2*p3*p4);
a1 = -real ( p2*p3*p4 + p1*p3*p4 + p1*p2*p4 + p1*p2*p3);
a2 = real (p1*p2 + p1*p3 + p1*p4 + p2*p3 + p2*p4 + p3*p4);
a3 = -real (p1 + p2 + p3 + p4);


w0 = sqrt(a2-a0*a3/a1);
w1 = sqrt(a0)/w0;
disp(['f0 = ', num2str(w0/2/pi)]);
disp(['f1 = ', num2str(w1/2/pi)]);
disp(['zeta1 = ', num2str(a1 / w0^2 / (2*w1))]);
disp(['mu = ', num2str(a3/a1*w0^2 - 1)]);


