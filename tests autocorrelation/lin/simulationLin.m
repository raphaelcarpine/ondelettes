clear all
close all


%% système

w0 = [10, 10.05, 20]*2*pi;
zeta = [0.01, 0.015, 0.01];
w0 = [10]*2*pi;
zeta = [0.01];

C = [1, -2, 1];
C = 1;

nbDDL = 3;


%% système

m = [1, 0.01, 2];
k = (10*2*pi)^2 * [1, 0.02, 0.02, 1];
c = 1 * [0, 0, 0, 1];
c = 0.1*[1 1 1 1];

M = diag(m);
K = zeros(3);
for ik = 1:4
    if ik - 1 > 0
        K(ik-1, ik-1) = K(ik-1, ik-1) + k(ik);
    end
    if ik <= 3
        K(ik, ik) = K(ik, ik) + k(ik);
    end
    if ik - 1 > 0 && ik <= 3
        K(ik, ik-1) = K(ik, ik-1) - k(ik);
        K(ik-1, ik) = K(ik-1, ik) - k(ik);
    end
end
C = zeros(3);
for ic = 1:4
    if ic - 1 > 0
        C(ic-1, ic-1) = C(ic-1, ic-1) + c(ic);
    end
    if ic <= 3
        C(ic, ic) = C(ic, ic) + c(ic);
    end
    if ic - 1 > 0 && ic <= 3
        C(ic, ic-1) = C(ic, ic-1) - c(ic);
        C(ic-1, ic) = C(ic-1, ic) - c(ic);
    end
end

syst = systLin(M, K, C);

[poles, shapes] = syst.normalModes()
[poles, shapes] = syst.complexModes()
[Mb, Kb, Cb] = syst.modalDamping();
Kb^(-1/2)*Cb


%% excitation

T = 100;
dt = 0.005;
fe = 1/dt;

t = 0:dt:T;
nt = length(t);

tdirac = 0;
Edirac = 1; % energie du signal dirac (integrale carré)
Tsyst = mean(1 ./ (w0.*zeta));
Pnoise = 0.0001 * Edirac / Tsyst; % puissance du signal noise (integrale carré / T)

f = zeros(1, nt);
f(floor(tdirac/dt)+1) = sqrt(Edirac/dt);
%f = f + sqrt(Pnoise) * randn(1, nt);

if true % bruit parfait pour fourier
    f = exp(2i*pi*rand(1, nt)); % bruit
    f(1) = 1;
    for k = 2:nt-1
        f(end+2-k) = conj(f(k));
    end
    f = ifft(f);
end

%% reponse

% ondelette
Q = 20;
MaxRidges = 3;
MaxParallelRidges = 3;
fmin = 5;
fmax = 25;
MultiSignalMode = true;

ct = 3;
cf = 5;

% plots

if true % f plot
    fig = figure('Name', 'excitation');
    ax = axes(fig);
    pltf = plot(ax, t, f);
    xlabel(ax, 't');
    ylabel(ax, 'x');
    
%     WaveletMenu('WaveletPlot', pltf, 'fmin', fmin, 'fmax', fmax,...
%         'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
%         , 'CtEdgeEffects', ct);
end

% x = reponseSystNddl(t, f, nbDDL, w0, zeta, C);
ddlF = 2;
x = syst.response(f, dt, ddlF);

fig = figure('Name', ['zeta=', num2str(zeta)]);
ax = axes(fig);
plt = plot(ax, t, x);
xlabel(ax, 't');
ylabel(ax, 'x');

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
    , 'CtEdgeEffects', ct, 'MultiSignalMode', MultiSignalMode);


%%

X = x;


return

%% autocorrelation

Rx = xcorr(x, 'biased') / var(x);
Rx = Rx(ceil(length(Rx)/2):end);


%%% test
if false
    n = length(x);
    Rx = nan(1, n);
    for k = 0:n-1
        Rx(k+1) = mean(x(1:n-k) .* x(k+1:end));
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
end


fig = figure('Name', ['autocorrelation;zeta=', num2str(zeta)]);
ax = axes(fig);
plt = plot(ax, t, Rx);
xlabel(ax, 't');
ylabel(ax, 'Rx');

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax,...
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



