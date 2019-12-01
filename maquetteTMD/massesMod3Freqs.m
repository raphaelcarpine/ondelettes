load('maquetteTMD/donneesEleves/mData.mat');

Deltam = 1.7756;


%% calcul freqs

data = Groupe4_1_1;

t = data(:,1)';
a = data(:,2)';


plt = plot(t, a);


fmin = 0.5;
fmax = 8;
Nf = 300;
Q = 30;
maxR = 3;
maxParallel = 3;
MinModu = 5;


WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq',...
    Nf, 'Q', Q, 'MaxRidges', maxR, 'MaxParallelRidges', maxParallel, 'RidgeMinModu', MinModu);


%%

F1 = 0.975;
F2 = 5.65;
F3 = 6.64;


%%


data = Groupe4_avec_masse_3_1;

t = data(:,1)';
a = data(:,2)';


plt = plot(t, a);


fmin = 0.5;
fmax = 8;
Nf = 300;
Q = 30;
maxR = 3;
maxParallel = 3;
MinModu = 5;


WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq',...
    Nf, 'Q', Q, 'MaxRidges', maxR, 'MaxParallelRidges', maxParallel, 'RidgeMinModu', MinModu);


%%

f1 = 0.95;
f2 = 5.42;
f3 = 6.405;




%%


m0 = [10, 10, 100];


fun = @(m) getFreqs3ddl(F1, F2, F3, m(1), m(2), m(3), Deltam) - [f1, f2, f3];


m = lsqnonlin(fun, m0);

M = nan(100, 3);
df = 1/40;
df = 1/200;
for i = 1:100
    F12 = F1 + df * randn;
    F22 = F2 + df * randn;
    F32 = F3 + df * randn;
    f12 = f1 + df * randn;
    f22 = f2 + df * randn;
    f32 = f3 + df * randn;
    fun = @(m) getFreqs3ddl(F12, F22, F32, m(1), m(2), m(3), Deltam) - [f12, f22, f32];
    m2 = lsqnonlin(fun, m);
    M(i, :) = m2;
end

Dm1 = std(M(:,1));
Dm2 = std(M(:,2));
Dm3 = std(M(:,3));







disp(['m1 = ', num2str(m(1)), ' +- ', num2str(Dm1)]);
disp(['m2 = ', num2str(m(2)), ' +- ', num2str(Dm2)]);
disp(['m3 = ', num2str(m(3)), ' +- ', num2str(Dm3)]);





%% results

% m1 = 33.4422;
% m2 = 24.1550;
% m3 = 17.1477;




