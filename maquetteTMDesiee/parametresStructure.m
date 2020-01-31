load('maquetteTMDesiee/donnees/mData.mat');

affichage = true;


%% structure seule

mesure = structure;
mesure = transpose(mesure);
t = mesure(1, :);
a = mesure(2, :);

f0 = getFreq(t, a, [0.8, 1.2]);
T0 = t(end) - t(1);

if affichage % affichage
    fig  = figure;
    ax = axes(fig);
    plt = plot(ax, t, a);
    xlabel(ax, 't');
    ylabel(ax, 'a');
    
    fmin = 0.5;
    fmax = 10;
    NbFreq = 300;
    Q = 10;
    
    WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq, 'Q', Q);
end

%% structure masse ajoutée

mesure = structureMasseAjoutee;
mesure = transpose(mesure);
t = mesure(1, :);
a = mesure(2, :);

f1 = getFreq(t, a, [0.8, 1.2]);
T1 = t(end) - t(1);

if affichage % affichage
    fig  = figure;
    ax = axes(fig);
    plt = plot(ax, t, a);
    xlabel(ax, 't');
    ylabel(ax, 'a');
    
    fmin = 0.5;
    fmax = 10;
    NbFreq = 300;
    Q = 10;
    
    WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq, 'Q', Q);
end



%% calculs

alpha = 1/(2*sqrt(3)); %ecart type loi uniforme

Df0 = alpha/T0;
Df1 = alpha/T1;

mu = delta_m/((f0/f1).^2-1);

Dmu = mu * sqrt(( (2*Df0*f0/f1^2)^2 +  (2*Df1*f0/f1^2)^2   )/((f0/f1)^2-1)^2);





disp(['f0 = ', num2str(f0), '+-', num2str(Df0)]);
disp(['f0'' = ', num2str(f1), '+-', num2str(Df1)]);
disp(['mu = ', num2str(mu), '+-', num2str(Dmu)]);




























