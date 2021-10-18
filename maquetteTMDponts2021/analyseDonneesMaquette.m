getData();

[fctDefModale, fctDefModaleAnimation] = deformeeMaquette(captPos, captDir);

%%

figure;
plt = plot(T, X);
legend(channelNames);
selectLine();

fmin = 0.2;
fmax = 50;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'RemoveMean', true,...
    'RealShapePlot', fctDefModale, 'AnimatedShapePlot', fctDefModaleAnimation,...
    'FrequencyScale', 'log', 'FourierScale', 'log');