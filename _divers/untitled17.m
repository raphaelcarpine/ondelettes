MotherWavelet = 'cauchy';
Q = 6;
delta = 0;
N = 100;

%%
dt = 0.01;
t = -10:dt:10;

freqs = nan(N, length(t));

figure;
for k = 1:100
    x = sin(2*pi*5*t + 2*pi*rand) .* (t<0) + 10*sin(2*pi*5*t + delta) .* (t>=0);
    
%     figure;
%     plt = plot(t, x);
%     WaveletMenu('WaveletPlot', plt);
    
    ridge = RidgeExtract(t, x, Q, 2, 8, 300, 'MotherWavelet', MotherWavelet,...
        'NbMaxRidges', 1, 'MaxSlopeRidge', inf, 'NbMaxParallelRidges', inf);
    
    plot(ridge.time{1}, ridge.freq{1});
    hold on
    
    ni = find(t == ridge.time{1}(1));
    nf = find(t == ridge.time{1}(end));
    
    freqs(k, ni:nf) = ridge.freq{1};
end

figure;
plot(t, mean(freqs));