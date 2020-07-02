clear all
% close all


%% système

w0 = [1*2*pi, 2*2*pi];
zeta = [0.01, 0.03];

C = [1, 1;
    0.5, -2];
% C = [1, 3];

nbDDL = 2;


%% excitation

T = 100;
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
%     f = - sqrt(2) *  erfcinv(2*rand(1, nt));
    f = randn(1, nt) / sqrt(dt);
else
    error('');
end

%% reponse

x = reponseSystNddl(t, f, nbDDL, w0, zeta, C);
x = x + 0;


% ondelette
Q = 5;
MaxRidges = 1;
MaxParallelRidges = 1;
fmin = 0.5;
fmax = 3;
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



fig = figure('Name', ['zeta=', num2str(zeta)]);
ax = axes(fig);
plt = plot(ax, t, x);
xlabel(ax, 't');
ylabel(ax, 'x');

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
    , 'CtEdgeEffects', ct);

%% test

% 
% Rx = crossCorrelation(x);
% 
% Nsv = 2;
% [SVrx, SVvectrx] = svdCWT(t, Rx, fmin, fmax, NbFreq, Q, Nsv);
% 
% WvltPlot2(t, linspace(fmin, fmax, NbFreq), SVrx{1}, 'module', Q, ct, 'log10');
% WvltPlot2(t, linspace(fmin, fmax, NbFreq), SVrx{2}, 'module', Q, ct, 'log10');

return

%% autocorrelation

Rx = xcorr(x, 'biased');% / var(x);
Rx = Rx(ceil(length(Rx)/2):end);


%%% test
if false
    n = length(x);
    Rx = nan(size(x));
    for k = 0:n-1
        Rx(:, k+1) = sum(x(:, 1:n-k) .* x(:, k+1:end), 2);
%         Rx(:, k+1) = sum(x .* [x(:, k+1:end), x(:, 1:k)], 2);
    end
end
%%% test

if true
    Rf = xcorr(f, 'biased') / var(f);
    Rf = Rf(ceil(length(Rf)/2):end);
    fig = figure('Name', 'autocorrelation;f');
    ax = axes(fig);
    plt = plot(ax, t, Rf);
    xlabel(ax, 't');
    ylabel(ax, 'Rf');
    
    WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
        'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
        , 'CtEdgeEffects', ct);
end


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



