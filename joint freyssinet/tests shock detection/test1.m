%% data

t = linspace(-10, 25, 100000);
t0 = 0;
f1 = 6;
f2 = 100;
x = sin(2*pi*f1*(t-t0)) .* exp(-0.02*2*pi*f1*(t-t0)) .* (t >= t0);
x = x + sin(2*pi*f2*(t-t0)) .* exp(-0.02*2*pi*f2*(t-t0)) .* (t >= t0);
x = x + randn(size(x));


%% plot

fig = figure;
ax = axes(fig);
plt = plot(ax, t, x);


fmin = 1;
fmax = 10;
Q = 10;
FrequencyScale = 'lin';
WvltScale = 'lin';
MotherWavelet =  'morlet';

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MultiSignalMode', false,...
    'WvltScale', WvltScale, 'FrequencyScale', FrequencyScale, 'MotherWavelet', MotherWavelet);%,...