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

F1 = 6.64;
F2 = 5.65;


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

f1 = 6.405;
f2 = 5.42;




%%


m0 = [17, 24];


fun = @(m) getFreqs2ddl(F1, F2, m(1), m(2), Deltam) - [f1, f2];

lb = [0 0];
ub = 10*[100 100];

options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e4);


m = lsqnonlin(fun, m0, lb, ub, options);

M = nan(10, 2);
df = 1/40;

for i = 1:10
    F12 = F1 + df * randn;
    F22 = F2 + df * randn;
    f12 = f1 + df * randn;
    f22 = f2 + df * randn;
    fun = @(m) getFreqs2ddl(F12, F22, m(1), m(2), Deltam) - [f12, f22];
    m2 = lsqnonlin(fun, m, lb, ub, options);
    M(i, :) = m2;
end

Dm1 = std(M(:,1));
Dm2 = std(M(:,2));





disp(['m1 = ', num2str(m(1)), ' +- ', num2str(Dm1)]);
disp(['m2 = ', num2str(m(2)), ' +- ', num2str(Dm2)]);





%% results

% m1 = 33.4422;
% m2 = 24.1550;
% m3 = 17.1477;




