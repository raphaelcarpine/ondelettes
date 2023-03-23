fs = 20;
dt = 1/fs;
t0 = dt*(0:25*fs-1);
N0 = length(t0);
p = sin(2*pi*t0) .* (sin(2*pi*t0/30).*(t0<=15) + 0.5*sin(2*pi*(t0-14)/10).*(t0>=14 & t0 <= 19));
p = sin(2*pi*t0) .* (t0 <= 10);
figure;
plot(t0, p)


t = dt*(0:500*fs-1);
N = length(t);
c = zeros(size(t));
Npeaks = 10;
c(floor(N*rand(1, Npeaks))+1) = exp(randn(1, Npeaks));

x = conv(c, p);
x = x(N0:end);

figure;
plot(t, x);

%%

% Fp0 = zeros(size(t0));
% for k = 1:8
%     Fp0 = Fp0 + abs(fft(x((k-1)*N0+1:k*N0))).^2;
% end
% Fp0 = sqrt(Fp0 / 8);

Rx = xcorr(x);
Fp0 = sqrt(abs(fft(Rx(N-249:N+250))));

p0 = ifft(Fp0);

Fx = fft(x);
M = @(theta) sum(abs(ifft(Fx ./ (fft([ifft(exp(1i*theta).*Fp0), zeros(1, N-N0)])))));
% M = @(theta) - abs(ifft(Fx ./ (fft([ifft(exp(1i*theta).*Fp0), zeros(1, N-N0)])))).^2;

% M(zeros(1, N0))
% M(2*pi*rand(1, N0))

opt = optimoptions('lsqnonlin','MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e4);
theta = lsqnonlin(M, 2*pi*zeros(1, N0), [], [], opt);

p0 = ifft(exp(1i*theta).*Fp0);
figure;
plot(t0, p0);
c0 = ifft(fft(x) ./ fft([p0, zeros(1, N-N0)]));
% figure;
% plot(t, c0);