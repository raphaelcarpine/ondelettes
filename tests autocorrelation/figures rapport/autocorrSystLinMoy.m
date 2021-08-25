%% data

% time
N = 100000;
N0 = 2000;
Dt = 0.05;

% averaging
Nav = 100;

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

N2 = floor(N/Nav);

Rx = zeros(1, N);
for kav = 1:Nav
    deltaK = (kav-1)*N2;
    for k = 0:N2-1
        Rx(k+1) = Rx(k+1) + Dt * sum( x(deltaK+(1:N2-k)) .* x(deltaK+(1+k:N2)) );
    end
end
Rx = Rx / Nav;
Rx = [Rx(end:-1:2), Rx];

Rh2 = zeros(1, N2);
for kav = 1:Nav
    deltaK = (kav-1)*N2;
    for k = 0:N2-1
        Rh2(k+1) = Rh2(k+1) + Dt * sum( h(deltaK+(1:N2-k)) .* h(deltaK+(1+k:N2)) );
    end
end
Rh2 = Rh2 / Nav;
Rh2 = [Rh2(end:-1:2), Rh2];

n = -N+1:N-1;

Rx0 = Rx(floor(end/2));

figure;
plot(n, Rx / Rx0);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$R_x(n)$', 'interpreter', 'latex', 'FontSize', 13);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
YLimRx = max(abs(get(gca, 'YLim') )) * [-1 1];
ylim(YLimRx);


figure;
plot(n, Dt^2*sigma^2 * (N2-abs(n)).* Rh  / Rx0);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$(N-|n|)\Delta t^2\sigma^2 R_h(n)$', 'interpreter', 'latex', 'FontSize', 13);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ylim(YLimRx);


figure;
plot(n, (Rx - Dt^2*sigma^2 * (N2-abs(n)).* Rh)  / Rx0);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$R_x(n) - (N-|n|)\Delta t^2\sigma^2 R_h(n)$', 'interpreter', 'latex', 'FontSize', 13);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ylim(YLimRx);
xlim([-1500 1500]);
ylim([-0.2 0.2]);


%%

alpha = sqrt( 4/(-real(lambda(1))*N*Dt) + 2/abs(imag(lambda(1))*N*Dt) )





