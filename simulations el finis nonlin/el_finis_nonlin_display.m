function el_finis_nonlin_display(varargin)

dispStresses = 0;
dispFreqs = 0;
plotShapes = 0;
plotAnimation = 1;
plotPosition = 1;
plotSpeed = 0;
plotAcc = 1;

%%

filePath = getResultsFile(varargin{:});

load(filePath);

%% calcul, mise en forme matrices

% interpolation etc

getYtot = @(Y) [zeros(1, size(Y, 2)); Y; zeros(1, size(Y, 2))];

% DDL 1, 2, end-1, end
Ytot = getYtot(Y);
Vtot = getYtot(V);
Atot = getYtot(A);

% capteurs
pos_capteurs = L/2;
Ycapt = getYcapt2(Ytot, pos_capteurs, dx);
Vcapt = getYcapt2(Vtot, pos_capteurs, dx);
Acapt = getYcapt2(Atot, pos_capteurs, dx);

%% display

% freqs et def modales
plotDefModales(dispStresses, dispFreqs, plotShapes, L, E, J, mu, N, Mr, Cr, Kr, fleche_pont, sigma0, sigma_min, sigma_max, delta_sigma_1t);

% animation
if plotAnimation
    time_coeff = 1.;
    movingPlotVehicles(Ytot, T, L, [t_vehicles_left, t_PL_left], [t_vehicles_right, t_PL_right],...
        [m_vehicles_left, m_PL_left], [m_vehicles_right, m_PL_right],...
        [c_vehicles_left, c_PL_left], [c_vehicles_right, c_PL_right],...
        pos_capteurs, x_nonlin, nonlin_reached, time_coeff);
end

% plots
fmin = 1;
fmax = 3;
Q = 2.5;
MaxRidges = 1;
MultiSignal = true;
[fctDefModale, fctDefModaleAnimation] = defModalePontLin(L, pos_capteurs);

if plotPosition
    figure;
    plt = plot(T, Ycapt);
    xlabel('Time [s]');
    ylabel('Displacement [m]');
    
    WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', MultiSignal,...
        'MaxParallelRidges', 1, 'MaxSlopeRidge', inf,...
        'MaxRidges', MaxRidges, 'RealShapePlot', fctDefModale, 'AnimatedShapePlot', fctDefModaleAnimation,...
        'SignalUnit', 'm', 'SquaredSignalUnit', 'm²');
end

if plotSpeed
    figure;
    plt = plot(T, Vcapt);
    xlabel('Time [s]');
    ylabel('Speed [m/s]');
    
    WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', MultiSignal,...
        'MaxParallelRidges', 1, 'MaxSlopeRidge', inf,...
        'MaxRidges', MaxRidges, 'RealShapePlot', fctDefModale, 'AnimatedShapePlot', fctDefModaleAnimation,...
        'SignalUnit', 'm/s', 'SquaredSignalUnit', 'm²/s²');
end

if plotAcc
    figure;
    plt = plot(T, Acapt);
    xlabel('Time [s]');
    ylabel('Acceleration [m/s²]');
    
    WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', MultiSignal,...
        'MaxParallelRidges', 1, 'MaxSlopeRidge', inf,...
        'MaxRidges', MaxRidges, 'RealShapePlot', fctDefModale, 'AnimatedShapePlot', fctDefModaleAnimation);
end


end