% test autocorr suite diracs

T = 1000;
dt = 0.01;

N = 1; %bd chocks
t0 = T*rand(1, N); t0 = 0;
nt0 = floor(t0/dt) + 1;
ampl0 = rand(1, N);
ampl0 = ones(1, N);
ampl0 = ampl0/dt;
ampl0 = sqrt(10)*ampl0;

t = 0:dt:T;
x = zeros(size(t));
x(nt0) = ampl0;

x = x + 2*randn(size(t));
x = rand(size(t));

% x = sin(2*pi*5*t).*(t>=100 & t<=200) + 2*sin(2*pi*5*t).*(t>=500 & t<=600);


%% plot

figure;
plt = plot(t, x);

WaveletMenu('WaveletPlot', plt);