load('maquetteTMD/donnees2/mData.mat');



%% frequence propre TMD

mesures = {TMDnonamorti11, TMDnonamorti12, TMDnonamorti13};
f11 = nan(1, 3);
for k = 1:3
    mesure = transpose(mesures{k});
    t = mesure(1,:);
    a = mesure(2,:);
    f11(k) = getFreq(t, a, [6.4, 7], 0.01);
end
Df11 = std(f11)/sqrt(3);
f11 = mean(f11);

mesures = {TMDnonamorti21, TMDnonamorti22};
f12 = nan(1, 2);
for k = 1:2
    mesure = transpose(mesures{k});
    t = mesure(1,:);
    a = mesure(2,:);
    f12(k) = getFreq(t, a, [6.4, 7], 0.001);
end
Df12 = std(f12)/sqrt(2);
f12 = mean(f12);

disp(['EPD : f1 = ', num2str(f11), '+-', num2str(Df11)]);
disp(['PBD : f1 = ', num2str(f12), '+-', num2str(Df12)]);

f1 = [f11, f12];


%% amortissement TMD

%t0 = nan(2, 5);
t0 = [1.2290    0.3850    0.5320    0.3980       NaN
    0.8830    0.4400    0.6430    0.6000    0.6030];

datas = {{TMD11, TMD12, TMD13, TMD14}, {TMD21, TMD22, TMD23, TMD24, TMD25}};

zeta1 = {nan(1, 4), nan(1, 5)};

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
        
        
        fmin = 5;
        fmax = 9;
        Nf = 300;
        Q = 10;
        maxR = 1;
        maxParallel = 1;
        
        
        WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax,...
            'NbFreq', Nf, 'Q', Q, 'MaxRidges', maxR, 'MaxParallelRidges', maxParallel);
        
        % regression
        A1lambda1 = RegressionMenu('FigureName', 'lambda2', 'Equation', 'A2*exp(-lambda2*x)',...
            'Param', 'A2 lambda2', 'Param0', [1 1], 'Fit', 'log(y)');
        lambda1= A1lambda1(2);
        zeta1{reglage}(kmesure) = lambda1/(2*pi*f1(reglage));
        
        delete(fig);
    end
end

% moyenne et ecart type
Dzeta11 = std(zeta1{1})/sqrt(length(zeta1{1}));
zeta11 = mean(zeta1{1});

Dzeta12 = std(zeta1{2})/sqrt(length(zeta1{2}));
zeta12 = mean(zeta1{2});

% Dzeta1 = [Dzeta11, Dzeta12];
% zeta1 = [zeta11, zeta12];

disp(['EPD : zeta1 = ', num2str(zeta11), '+-', num2str(Dzeta11)]);
disp(['PBD : zeta1 = ', num2str(zeta12), '+-', num2str(Dzeta12)]);




























