%% data

% time
N = 100000;
N0 = 2000;
Dt = 0.05;

% system
lambda = -0.1 + 2i*pi;
A = -1i;
lambda = [lambda, conj(lambda)];
A = [A, conj(A)];
h = sum(A.' .* exp(Dt * lambda.' * (0:N+N0-1)), 1);
% figure;
% plot(h);
B = Dt * sum(A.' * conj(A) ./ ( 1 - exp(Dt*lambda).' * conj(exp(Dt*lambda))), 1);
Rh = sum(B.' .* exp(Dt * lambda.' * (0:N-1)), 1);
Rh = [Rh(end:-1:2), Rh];
% figure;
% plot(Rh);

% noise
sigma = 1;

%% response

w = sigma*randn(1, N + N0);

x = nan(1, N);
for k = 1:N % convolution
    x(k) = Dt * sum( w(1:N0+k) .* h(N0+k:-1:1) );
end
% figure;
% plot(x);

%% autocorr

Rx = nan(1, N);
for k = 0:N-1
    Rx(k+1) = Dt * sum( x(1:N-k) .* x(1+k:N) );
end
Rx = [Rx(end:-1:2), Rx];

Rh2 = nan(1, N);
for k = 0:N-1
    Rh2(k+1) = Dt * sum( h(1:N-k) .* h(1+k:N) );
end
Rh2 = [Rh2(end:-1:2), Rh2];

alpha = sqrt( 4/(-real(lambda(1))*N*Dt) + 2/abs(imag(lambda(1))*N*Dt))
A0 = Dt^2*sigma^2*N*Rh((end-1)/2+1);


n = -N+1:N-1;

I = n >= -5000 & n <= 5000;
n = n(I);
Rx = Rx(I);
Rh = Rh(I);

Nmax = 3000;

figure;
plt1 = plot(n, Rx);
yline(alpha*A0, '--r');
yline(-alpha*A0, '--r');
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 15);
ylabel('$R_{\tilde x}$', 'interpreter', 'latex', 'FontSize', 15);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
YLimRx = max(abs(get(gca, 'YLim') )) * [-1 1];
ylim(YLimRx);
xlim([-Nmax, Nmax]);


figure;
plt2 = plot(n, Dt^2*sigma^2 * (N-abs(n)).* Rh);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 15);
ylabel('$(N-|n|)\Delta t^2\sigma^2 R_h$', 'interpreter', 'latex', 'FontSize', 15);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ylim(YLimRx);
xlim([-Nmax, Nmax]);


figure;
plt3 = plot(n, Rx - Dt^2*sigma^2 * (N-abs(n)).* Rh);
yline(alpha*A0, '--r');
yline(-alpha*A0, '--r');
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 15);
ylabel('$R_{\tilde x} - (N-|n|)\Delta t^2\sigma^2 R_h$', 'interpreter', 'latex', 'FontSize', 15);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ylim(YLimRx);
xlim([-Nmax, Nmax]);


%%


%% CWT

WaveletMenu('WaveletPlot', plt1, 'fmin', 0.03, 'fmax', 0.07);
WaveletMenu('WaveletPlot', plt2, 'fmin', 0.03, 'fmax', 0.07);
WaveletMenu('WaveletPlot', plt3, 'fmin', 0.03, 'fmax', 0.07);



