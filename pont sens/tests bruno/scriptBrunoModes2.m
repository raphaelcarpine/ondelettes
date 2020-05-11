% choix precision
ct = 3;
cf = 5;


% construction du signal
f1 = 10;
f2 = 12;
zeta1 = 0.01;
zeta2 = 0.01;
shape1 = [1; 3; -1; 2];
shape2 = [1; -2; -1; 1+0.1i];
T = 20;

t = linspace(0, T, 10000);
X = shape1 * exp(sqrt(1-zeta1^2)*2i*pi*f1*t - zeta1*2*pi*f1*t + 0.2i*pi*sin(t));
X = X + shape2 * exp(sqrt(1-zeta2^2)*2i*pi*f2*t - zeta2*2*pi*f2*t);
X = real(X);

% choix Q
Dt1 = 1/(zeta1*2*pi*f1);
Df1 = f2-f1;
[Qmin, Qmax, Qz] = getBoundsQ(f1, Df1, Dt1, T, ct, cf);
disp('premier mode :');
disp(['Qmin = ', num2str(Qmin), ' ; Qmax = ', num2str(Qmax), ' ; Qz = ', num2str(Qz)]);

Dt2 = 1/(zeta2*2*pi*f2);
Df2 = f2-f1;
[Qmin, Qmax, Qz] = getBoundsQ(f2, Df2, Dt2, T, ct, cf);
disp('deuxième mode :');
disp(['Qmin = ', num2str(Qmin), ' ; Qmax = ', num2str(Qmax), ' ; Qz = ', num2str(Qz)]);

Q = 40;
disp(' ');
disp(['Q = ', num2str(Q)]);
disp(' ');

% shape plot
bridgeDim = [17, 4];
channelPos = [3, 1.5; 9, 2; 12, 2; 16, 2.5];
shapePlotBridge = @(shape, figTitle) shapePlotPlate(bridgeDim, channelPos, shape, figTitle);

% CWT
NbMaxRidges = 2;
NbMaxParallelRidges = inf;
fmin = 8;
fmax = 14;
nbFreqs = 300;

figure;
plotX = plot(t, X);

WaveletMenu('WaveletPlot', plotX, 'Q', Q, 'fmin', fmin, 'fmax', fmax, 'NbFreq', nbFreqs, 'MultiSignalMode', true,...
    'MaxRidges', NbMaxRidges, 'MaxParallelRidges', NbMaxParallelRidges, 'RealShapePlot', shapePlotBridge);
