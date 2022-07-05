dataFile = 'vide_p1.csv';

getData();

% [fctDefModale, fctDefModaleAnimation] = deformeeMaquette(captPos, captDir);

%%
close all

figure;
plt = plot(T, X);
legend(channelNames);
selectLine();

fmin = 0.1;
fmax = 1.5;
Q = 3;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'RemoveMean', true);%,...
%     'RealShapePlot', fctDefModale, 'AnimatedShapePlot', fctDefModaleAnimation,...
%     'FrequencyScale', 'lin', 'FourierScale', 'lin');