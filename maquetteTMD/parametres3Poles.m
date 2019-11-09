load('maquetteTMD/donnees2/mData.mat');


%% calcul poles

t0 = nan(2, 6);
t0 = [0.7680    1.6100    0.8420    0.7550    0.9550       NaN
    1.8300    1.0780    0.8800    0.9800    0.8500    1.0600];

datas = {{StructureTMD11, StructureTMD12, StructureTMD13, StructureTMD14, StructureTMD15},...
    {StructureTMD21, StructureTMD22, StructureTMD23, StructureTMD24, StructureTMD25, StructureTMD26}};

f0 = {nan(1, 5), nan(1, 6)};
f02 = {nan(1, 5), nan(1, 6)};
f1 = {nan(1, 5), nan(1, 6)};
zeta1 = {nan(1, 5), nan(1, 6)};
mu = {nan(1, 5), nan(1, 6)};
mu2 = {nan(1, 5), nan(1, 6)};

reglage = 1;

data = datas{reglage};

for kmesure = 1:length(data)
    mesure = data{kmesure};
    
    mesure = transpose(mesure);
    t = mesure(1,:);
    t = t-t(1);
    a = mesure(2,:);
    
    % affichage pour déterminer t0
    if isnan(t0(reglage, kmesure))
        fig = figure;
        ax = axes(fig);
        plot(ax, t, a);
        
        t0(reglage, kmesure) = input('t0 = ');
        
        delete(fig);
    end
    
    a = a(t>=t0(reglage, kmesure));
    t = t(t>=t0(reglage, kmesure));
    
    % wavelet
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
    
    % regression
    fa = RegressionMenu('FigureName', 'fa', 'Equation', 'fa', 'Param', 'fa', 'Param0', 6);
    fb = RegressionMenu('FigureName', 'fb', 'Equation', 'fb', 'Param', 'fb', 'Param0', 6);
    fc = RegressionMenu('FigureName', 'fc', 'Equation', 'fc', 'Param', 'fc', 'Param0', 6);
    
    Aala = RegressionMenu('FigureName', 'la', 'Equation', 'Aa*exp(-la*x)',...
        'Param', 'Aa la', 'Param0', [1, 1], 'Fit', 'log(y)');
    Ablb = RegressionMenu('FigureName', 'lb', 'Equation', 'Ab*exp(-lb*x)',...
        'Param', 'Ab lb', 'Param0', Aala, 'Fit', 'log(y)');
    Aclc = RegressionMenu('FigureName', 'lc', 'Equation', 'Ac*exp(-lc*x)',...
        'Param', 'Ac lc', 'Param0', Ablb, 'Fit', 'log(y)');
    la = Aala(2);
    lb = Ablb(2);
    lc = Aclc(2);
    
    % parametres systeme
    p1 = -la + 1i*2*pi*fa;
    p2 = -la - 1i*2*pi*fa;
    p3 = -lb + 1i*2*pi*fb;
    p4 = -lb - 1i*2*pi*fb;
    p5 = -lc + 1i*2*pi*fc;
    p6 = -lc - 1i*2*pi*fc;
    
    poles = sort([p1, p2, p3, p4, p5, p6]);
    
    
    S = @(param) abs(poles - getPoles3ddl(param(1), param(2), param(3), param(4), param(5), param(6)));
    
    if kmesure == 1
        param0 = [6.65, 5.7, 0.014, 0.005, 6.55, 0.07];
    else
        param0 = [f0{reglage}(kmesure-1), f02{reglage}(kmesure-1), mu{reglage}(kmesure-1),...
            mu2{reglage}(kmesure-1), f1{reglage}(kmesure-1), zeta1{reglage}(kmesure-1)];
    end
    
    
    optionsReg = optimoptions(@lsqnonlin, 'MaxIterations', 1e4,...
        'StepTolerance', 1e-6, 'MaxFunctionEvaluations', inf, 'FunctionTolerance', 0);
    lbound = zeros(1, 6);
    ubound = inf*ones(1, 6);
    
    param = lsqnonlin(S, param0, lbound, ubound, optionsReg);
    
    
    
    f0{reglage}(kmesure) = param(1);
    f02{reglage}(kmesure) = param(2);
    mu{reglage}(kmesure) = param(3);
    mu2{reglage}(kmesure) = param(4);
    f1{reglage}(kmesure) = param(5);
    zeta1{reglage}(kmesure) = param(6);
    
    % fermeture de la fig
    try
        close(fig);
    catch
    end
end

% moyenne et ecart type
Df0 = std(f0{reglage})/sqrt(length(f0{reglage}));
f0 = mean(f0{reglage});

Df1 = std(f1{reglage})/sqrt(length(f1{reglage}));
f1 = mean(f1{reglage});

Df02 = std(f02{reglage})/sqrt(length(f02{reglage}));
f02 = mean(f02{reglage});

Dzeta1 = std(zeta1{reglage})/sqrt(length(zeta1{reglage}));
zeta1 = mean(zeta1{reglage});

Dmu = std(mu{reglage})/sqrt(length(mu{reglage}));
mu = mean(mu{reglage});

Dmu2 = std(mu2{reglage})/sqrt(length(mu2{reglage}));
mu2 = mean(mu2{reglage});


disp(['f0 = ', num2str(f0), '+-', num2str(Df0)]);
disp(['f02 = ', num2str(f02), '+-', num2str(Df02)]);
disp(['f1 = ', num2str(f1), '+-', num2str(Df1)]);
disp(['zeta1 = ', num2str(100*zeta1), '+-', num2str(100*Dzeta1), ' %']);
disp(['mu = ', num2str(mu), '+-', num2str(Dmu)]);
disp(['mu2 = ', num2str(mu2), '+-', num2str(Dmu2)]);

%% results







