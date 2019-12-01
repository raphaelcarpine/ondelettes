%% chargement des données

load('maquetteTMD/donneesTP/mData.mat');

data = TMD_structure_accelero_8;
% data = nan; % à compléter


mesure = transpose(data);
t = mesure(1,:);
t = t-t(1);
a = mesure(3,:);

%% affichage pour déterminer t0

fig = figure;
ax = axes(fig);
plot(ax, t, a);

t0 = input('t0 = ');

delete(fig);

a = a(t>=t0);
t = t(t>=t0);

%% transformée en ondelettes

fig = figure;
ax = axes(fig);
plt = plot(ax, t, a);


fmin = 5;
fmax = 8;
Nf = 300;
Q = 30;
maxR = 3;
maxParallel = 3;
MinModu = 5;


WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq',...
    Nf, 'Q', Q, 'MaxRidges', maxR, 'MaxParallelRidges', maxParallel, 'RidgeMinModu', MinModu);

%% régression pour déterminer les pôles 

fa = RegressionMenu('FigureName', 'fa', 'Equation', 'fa', 'Param', 'fa', 'Param0', 6);
fb = RegressionMenu('FigureName', 'fb', 'Equation', 'fb', 'Param', 'fb', 'Param0', 6);
fc = RegressionMenu('FigureName', 'fc', 'Equation', 'fc', 'Param', 'fc', 'Param0', 6);

Aalambdaa = RegressionMenu('FigureName', 'lambdaa', 'Equation', 'Aa*exp(-lambdaa*x)',...
    'Param', 'Aa lambdaa', 'Param0', [1, 1], 'Fit', 'log(y)');
Ablambdab = RegressionMenu('FigureName', 'lambdab', 'Equation', 'Ab*exp(-lambdab*x)',...
    'Param', 'Ab lambdab', 'Param0', Aalambdaa, 'Fit', 'log(y)');
Aclambdac = RegressionMenu('FigureName', 'lambdac', 'Equation', 'Ac*exp(-lambdac*x)',...
    'Param', 'Ac lambdac', 'Param0', Ablambdab, 'Fit', 'log(y)');
lambdaa = Aalambdaa(2);
lambdab = Ablambdab(2);
lambdac = Aclambdac(2);

%%

pa = nan; % à compléter
pb = nan; % à compléter
pc = nan; % à compléter

pa = -lambdaa + 1i*2*pi*fa;
pb = -lambdab + 1i*2*pi*fb;
pc = -lambdac + 1i*2*pi*fc;

disp(['p_a = ', num2str(pa)]);
disp(['p_a* = ', num2str(conj(pa))]);
disp(['p_b = ', num2str(pb)]);
disp(['p_b* = ', num2str(conj(pb))]);
disp(['p_c = ', num2str(pc)]);
disp(['p_c* = ', num2str(conj(pc))]);


poles = sort([pa, conj(pa), pb, conj(pb), pc, conj(pc)]);




