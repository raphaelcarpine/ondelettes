f = 5;
% phi = 2*pi*rand();
dPhi = -0.1*2*pi;
% deltaTdPhi = 0.2 * 1/f * 0;
deltaTdPhi = dPhi / (2*pi*f);
phi = - dPhi/2;
Fs = 200;


if deltaTdPhi ~= 0
    H = @(t) dPhi*(t/deltaTdPhi + 1/2).*(t>=-deltaTdPhi/2) - dPhi*(t/deltaTdPhi - 1/2).*(t>=deltaTdPhi/2);
else
    H = @(t) dPhi*(t>=0);
end

t = -10:1/Fs:10;
x = cos(2*pi*f*t + phi + H(t));
% x = exp(1i*(2*pi*f*t + phi + H(t)));

if false % plot H(t)
    figure;
    plot(t, H(t));
end
if false % plot phase
    figure;
    plot(t, 2*pi*f*t + phi + H(t));
end

figure;
plt = plot(t, x);
ylim([-1.5, 1.5]);

%% wavelet

fmin = 0.01;
fmax = 100;
% fmin = 10;
% fmax = 50;
Q = 5;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MotherWavelet', 'morlet');
