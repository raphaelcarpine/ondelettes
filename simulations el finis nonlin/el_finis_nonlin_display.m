function el_finis_nonlin_display(varargin)

plotAnimation = 1;
plotPosition = 0;
plotSpeed = 0;
plotAcc = 1;

%%

filePath = getResultsFile(varargin{:});

load(filePath);

%% display

% animation
if plotAnimation
    time_coeff = 1.;
    movingPlotVehicles(Ytot, T, L, t_vehicles_left, t_vehicles_right, m_vehicles_left, m_vehicles_right,...
        c_vehicles_left, c_vehicles_right, pos_capteurs, x_nonlin, nonlin_reached, time_coeff);
end

% wavelet capteurs
Ycapt = getYcapt2(Ytot, pos_capteurs, dx);
Vcapt = getYcapt2(Vtot, pos_capteurs, dx);
Acapt = getYcapt2(Atot, pos_capteurs, dx);
%Ycapt = Ycapt + 1e-6 * randn(size(Ycapt));


fmin = 0;
fmax = 1000;
Q = 5;
MaxRidges = 1;
[fctDefModale, fctDefModaleAnimation] = defModalePontLin(L, pos_capteurs);

if plotPosition
    figure;
    plt = plot(T, Ycapt);
    
    WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', false,...
        'MaxRidges', MaxRidges, 'RealShapePlot', fctDefModale, 'AnimatedShapePlot', fctDefModaleAnimation);
end

if plotSpeed
    figure;
    plt = plot(T, Vcapt);
    WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', false,...
        'MaxRidges', MaxRidges, 'RealShapePlot', fctDefModale, 'AnimatedShapePlot', fctDefModaleAnimation);
end

if plotAcc
    figure;
    plt = plot(T, Acapt);
    WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', false,...
        'MaxRidges', MaxRidges, 'RealShapePlot', fctDefModale, 'AnimatedShapePlot', fctDefModaleAnimation);
end


end