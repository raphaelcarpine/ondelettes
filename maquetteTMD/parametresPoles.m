load('maquetteTMD/donnees2/mData.mat');


%% calcul poles

t0 = nan(2, 6);
t0 = [0.7680    1.6100    0.8420    0.7550    0.9550       NaN
    1.8300    1.0780    0.8800    0.9800    0.8500    1.0600];

datas = {{StructureTMD11, StructureTMD12, StructureTMD13, StructureTMD14, StructureTMD15},...
    {StructureTMD21, StructureTMD22, StructureTMD23, StructureTMD24, StructureTMD25, StructureTMD26}};

f0 = {nan(1, 5), nan(1, 6)};
f1 = {nan(1, 5), nan(1, 6)};
zeta1 = {nan(1, 5), nan(1, 6)};
mu = {nan(1, 5), nan(1, 6)};

for reglage = 1:2
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
        
        
        fmin = 6;
        fmax = 8;
        Nf = 300;
        Q = 30;
        maxR = 2;
        maxParallel = 2;
        MinModu = '10*noise';
        
        
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
        f0{reglage}(kmesure) = w0/2/pi;
        f1{reglage}(kmesure) = w1/2/pi;
        zeta1{reglage}(kmesure) = a1 / w0^2 / (2*w1);
        mu{reglage}(kmesure) = a3/a1*w0^2 - 1;
        
        % fermeture de la fig
        try
            close(fig);
        catch
        end
    end
end

% moyenne et ecart type
Df01 = std(f0{1})/sqrt(length(f0{1}));
f01 = mean(f0{1});
Df02 = std(f0{2})/sqrt(length(f0{2}));
f02 = mean(f0{2});

Df11 = std(f1{1})/sqrt(length(f1{1}));
f11 = mean(f1{1});
Df12 = std(f1{2})/sqrt(length(f1{2}));
f12 = mean(f1{2});

Dzeta11 = std(zeta1{1})/sqrt(length(zeta1{1}));
zeta11 = mean(zeta1{1});
Dzeta12 = std(zeta1{2})/sqrt(length(zeta1{2}));
zeta12 = mean(zeta1{2});

Dmu1 = std(mu{1})/sqrt(length(mu{1}));
mu1 = mean(mu{1});
Dmu2 = std(mu{2})/sqrt(length(mu{2}));
mu2 = mean(mu{2});


disp('EDP :');
disp(['f0 = ', num2str(f01), '+-', num2str(Df01)]);
disp(['f1 = ', num2str(f11), '+-', num2str(Df11)]);
disp(['zeta1 = ', num2str(100*zeta11), '+-', num2str(100*Dzeta11), ' %']);
disp(['mu = ', num2str(mu1), '+-', num2str(Dmu1)]);
disp('');
disp('PBD :');
disp(['f0 = ', num2str(f02), '+-', num2str(Df02)]);
disp(['f1 = ', num2str(f12), '+-', num2str(Df12)]);
disp(['zeta1 = ', num2str(100*zeta12), '+-', num2str(100*Dzeta12), ' %']);
disp(['mu = ', num2str(mu2), '+-', num2str(Dmu2)]);


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



%% test

mesure = mesure0;

%mesure = transpose(mesure);
t = mesure(1,:);
t = t-t(1);
a = mesure(2,:);

% affichage pour déterminer t0
fig = figure;
ax = axes(fig);
plot(ax, t, a);

t0 = input('t0 = ');

delete(fig);

a = a(t>=t0);
t = t(t>=t0);

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
MinModu = 0;


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





