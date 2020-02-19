clear all
close all


%% système

w0 = 2*pi;
zeta = 0.02;

w02 = 3*pi;
zeta2 = 0.01;
C1 = 1;
C2 = 1;

nbDDL = 1;


%% excitation

T = 1000;
dt = 0.01;
fe = 1/dt;

t = 0:dt:T;
nt = length(t);

excitation = 'gaussien'; % 'bruit' 'dirac' 'gaussien'

if isequal(excitation, 'bruit')
    f = exp(2i*pi*rand(1, nt)); % bruit
    f(1) = 1;
    for k = 2:nt-1
        f(end+2-k) = conj(f(k));
    end
    f = ifft(f);
elseif isequal(excitation, 'dirac')
    f = ones(1, nt); % dirac
    f = ifft(f);
elseif isequal(excitation, 'gaussien')
%     f = normrnd(0, 1, 1, nt);
    f = rand(1, nt);
else
    error('');
end

%% reponse

% ondelette
Q = 30;
MaxRidges = 1;
MaxParallelRidges = 1;
fmin = 0.1;
fmax = 2;
NbFreq = 300;

ct = 3;
cf = 5;

% plots

if true % f plot
    fig = figure('Name', 'excitation');
    ax = axes(fig);
    pltf = plot(ax, t, f);
    xlabel(ax, 't');
    ylabel(ax, 'x');
    
    
    WaveletMenu('WaveletPlot', pltf, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
        'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
        , 'CtEdgeEffects', ct);
end

if nbDDL == 1
    x = reponseSyst1ddl(t, f, w0, zeta);
elseif nbDDL == 2
    x = reponseSyst2ddl(t, f, w0, w02, zeta, zeta2, C1, C2);
end

fig = figure('Name', ['zeta=', num2str(zeta)]);
ax = axes(fig);
plt = plot(ax, t, x);
xlabel(ax, 't');
ylabel(ax, 'x');

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
    , 'CtEdgeEffects', ct);


%% autocorrelation

% x = f;

Rx = xcov(x) / var(x);
Rx = Rx(ceil(length(Rx)/2):end);


%%% test
n = length(x);
Rx = nan(1, n);
for k = 0:n-1
    Rx(k+1) = mean(x(1:n-k) .* x(k+1:end));
end
%%% test


fig = figure('Name', ['autocorrelation;zeta=', num2str(zeta)]);
ax = axes(fig);
plt = plot(ax, t, Rx);
xlabel(ax, 't');
ylabel(ax, 'Rx');

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
    , 'CtEdgeEffects', ct);



return

%% stabilisation

freqs = (0:(length(x)-1)) / (t(end)-t(1));
freqs = freqs(1:floor(length(freqs)/2));
fftX = fft(x);
fftX = fftX(1:floor(length(fftX)/2));

figure;
modalsd(fftX.', freqs, fe, 'MaxModes', 20);



