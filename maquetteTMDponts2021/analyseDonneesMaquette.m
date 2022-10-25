dataFile = 'SensorConnectDataT100.csv';

getData();

[fctDefModale, fctDefModaleAnimation] = deformeeMaquette(captPos, captDir);

%%

figure;
plt = plot(T, X);
legend(channelNames);
selectLine();

fmin = 0.;
fmax = 2;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'RemoveMean', true,...
    'RealShapePlot', fctDefModale, 'AnimatedShapePlot', fctDefModaleAnimation,...
    'FrequencyScale', 'lin', 'FourierScale', 'lin');