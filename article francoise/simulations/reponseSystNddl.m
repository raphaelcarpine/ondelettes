function x = reponseSystNddl(t, f, N, w0, zeta, C)
%REPONSESYST2DDL Système linéaire de degré N avec amortissement
%proportionnel
%   Ck = (phik)^2/mk



%% equations

% z1'' + 2*zeta1*w01*z1' + w01^2*z1 = f * phi1/m1
% z2'' + 2*zeta2*w02*z2' + w02^2*z2 = f * phi2/m2
% ...
% x = phi1*z1 + phi2*z2 + ...

    function h = H(freq)
        h = 0;
        for k = 1:N
            h = h + C(k) ./ (-(2*pi*freq).^2 + 2i*zeta(k)*w0(k)*2*pi*freq + w0(k)^2);
        end
    end

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

