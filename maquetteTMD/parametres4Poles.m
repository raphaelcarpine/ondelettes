load('maquetteTMD/donnees2/mData.mat');


%% calcul poles

t0 = nan(2, 6);
t0 = [0.7680    1.6100    0.8420    0.7550    0.9550       NaN
    1.8300    1.0780    0.8800    0.9800    0.8500    1.0600];

datas = {{StructureTMD11, StructureTMD12, StructureTMD13, StructureTMD14, StructureTMD15},...
    {StructureTMD21, StructureTMD22, StructureTMD23, StructureTMD24, StructureTMD25, StructureTMD26}};

f0 = {nan(1, 5), nan(1, 6)};
f02 = {nan(1, 5), nan(1, 6)};
f03 = {nan(1, 5), nan(1, 6)};
f1 = {nan(1, 5), nan(1, 6)};
zeta1 = {nan(1, 5), nan(1, 6)};
mu = {nan(1, 5), nan(1, 6)};
mu2 = {nan(1, 5), nan(1, 6)};
mu3 = {nan(1, 5), nan(1, 6)};

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
    fd = RegressionMenu('FigureName', 'fd', 'Equation', 'fd', 'Param', 'fd', 'Param0', 6);
    
    Aala = RegressionMenu('FigureName', 'la', 'Equation', 'Aa*exp(-la*x)',...
        'Param', 'Aa la', 'Param0', [1, 1], 'Fit', 'log(y)');
    Ablb = RegressionMenu('FigureName', 'lb', 'Equation', 'Ab*exp(-lb*x)',...
        'Param', 'Ab lb', 'Param0', Aala, 'Fit', 'log(y)');
    Aclc = RegressionMenu('FigureName', 'lc', 'Equation', 'Ac*exp(-lc*x)',...
        'Param', 'Ac lc', 'Param0', Ablb, 'Fit', 'log(y)');
    Adld = RegressionMenu('FigureName', 'ld', 'Equation', 'Ad*exp(-ld*x)',...
        'Param', 'Ad ld', 'Param0', Aclc, 'Fit', 'log(y)');
    la = Aala(2);
    lb = Ablb(2);
    lc = Aclc(2);
    ld = Adld(2);
    
    % parametres systeme
    p1 = -la + 1i*2*pi*fa;
    p2 = -la - 1i*2*pi*fa;
    p3 = -lb + 1i*2*pi*fb;
    p4 = -lb - 1i*2*pi*fb;
    p5 = -lc + 1i*2*pi*fc;
    p6 = -lc - 1i*2*pi*fc;
    p7 = -ld + 1i*2*pi*fd;
    p8 = -ld - 1i*2*pi*fd;
    
    poles = sort([p1, p2, p3, p4, p5, p6, p7, p8]);
    
    
    S = @(param) abs(poles - getPoles4ddl(param(1), param(2), param(3), param(4),...
        param(5), param(6), param(7), param(8)));
    
    if kmesure == 1
        param0 = [6.65, 5.7, 1, 0.014, 0.01, 0.01, 6.55, 0.07];
    else
        param0 = param;
    end
    
    
    optionsReg = optimoptions(@lsqnonlin, 'MaxIterations', 1e5,...
        'StepTolerance', 1e-6, 'MaxFunctionEvaluations', inf, 'FunctionTolerance', 0);
    lbound = zeros(1, 8);
    ubound = inf*ones(1, 8);
    
    param = lsqnonlin(S, param0, lbound, ubound, optionsReg);
    
    
    
    f0{reglage}(kmesure) = param(1);
    f02{reglage}(kmesure) = param(2);
    f03{reglage}(kmesure) = param(3);
    mu{reglage}(kmesure) = param(4);
    mu2{reglage}(kmesure) = param(5);
    mu3{reglage}(kmesure) = param(6);
    f1{reglage}(kmesure) = param(7);
    zeta1{reglage}(kmesure) = param(8);
    
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

Df03 = std(f03{reglage})/sqrt(length(f03{reglage}));
f03 = mean(f03{reglage});

Dzeta1 = std(zeta1{reglage})/sqrt(length(zeta1{reglage}));
zeta1 = mean(zeta1{reglage});

Dmu = std(mu{reglage})/sqrt(length(mu{reglage}));
mu = mean(mu{reglage});

Dmu2 = std(mu2{reglage})/sqrt(length(mu2{reglage}));
mu2 = mean(mu2{reglage});

Dmu3 = std(mu3{reglage})/sqrt(length(mu3{reglage}));
mu3 = mean(mu3{reglage});


disp(['f0 = ', num2str(f0), '+-', num2str(Df0)]);
disp(['f02 = ', num2str(f02), '+-', num2str(Df02)]);
disp(['f03 = ', num2str(f03), '+-', num2str(Df03)]);
disp(['f1 = ', num2str(f1), '+-', num2str(Df1)]);
disp(['zeta1 = ', num2str(100*zeta1), '+-', num2str(100*Dzeta1), ' %']);
disp(['mu = ', num2str(mu), '+-', num2str(Dmu)]);
disp(['mu2 = ', num2str(mu2), '+-', num2str(Dmu2)]);
disp(['mu3 = ', num2str(mu3), '+-', num2str(Dmu3)]);


%% results

% EDP :
% f0 = 6.6342+-0.0012314
% f1 = 6.8686+-0.004029
% zeta1 = 6.2919+-0.045224 %
% mu = 0.011446+-7.1325e-05
% PBD :
% f0 = 6.6525+-0.0010866
% f1 = 6.802+-0.0027288
% zeta1 = 7.4354+-0.029154 %
% mu = 0.0098506+-2.4902e-05





