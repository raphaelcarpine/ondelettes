function x = reponseSyst1ddl(t, f, w0, zeta)
%REPONSESYST1DDL Summary of this function goes here
%   Detailed explanation goes here

%% equations
"x'' + 2*zeta*w0*x' + w0^2*x = f";

H = @(freq) 1 ./ (-(2*pi*freq).^2 + 2i*zeta*w0*2*pi*freq + w0^2);

%% resolution
dt = mean(diff(t));
if max(abs(diff(t)-dt)/dt) > 1e-2
    error('pas de temps non constant');
end

TFf = fft(f);
freqs = (0:(length(t)-1)) / (t(end)-t(1));


TFx = TFf .* H(freqs);

x = ifft(TFx);
x = real(x);

end

