%% figures difference R_{\tilde x} R_{h*\tilde w}

%% data

% time
N = 2000;
N1 = 1000;
N0 = 2000;
Dt = 0.05;

% system
lambda = -0.1 + 2i*pi;
A = -1i;
lambda = [lambda, conj(lambda)];
A = [A, conj(A)];
h = sum(A.' .* exp(Dt * lambda.' * (0:N+2*N1+N0-1)), 1);
% figure;
% plot(h);

% noise
sigma = 1;

%% x1

w = sigma*randn(1, N + 2*N1 + N0);

x = nan(1, N+2*N1);
for k = 1:N+2*N1
    x(k) = Dt * sum( w(1:N0+k) .* h(N0+k:-1:1) );
end
x1 = x .* [zeros(1, N1), ones(1, N), zeros(1, N1)];

n = -N1:N+N1-1;

figure;
plot(n, w(N0+1:end));
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$w$', 'interpreter', 'latex', 'FontSize', 13);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
YLimW = max(abs(get(gca, 'YLim') )) * [-1 1];
ylim(YLimW);

figure;
plot(n, x, ':');
hold on
set(gca,'ColorOrderIndex', get(gca,'ColorOrderIndex')-1)
plot(n, x1);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$\tilde x$', 'interpreter', 'latex', 'FontSize', 13);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
YLimX = max(abs(get(gca, 'YLim') )) * [-1 1];
ylim(YLimX);

%% x2

w2 = w .* [zeros(1, N0+N1), ones(1, N), zeros(1, N1)];

x2 = nan(1, N+2*N1);
for k = 1:N+2*N1
    x2(k) = Dt * sum( w2(1:N0+k) .* h(N0+k:-1:1) );
end


figure;
plot(n, w2(N0+1:end));
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$\tilde w$', 'interpreter', 'latex', 'FontSize', 13);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ylim(YLimW);

figure;
plot(n, x2);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$h*\tilde w$', 'interpreter', 'latex', 'FontSize', 13);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ylim(YLimX);

%% x2 - x1

figure;
plot(n, x2-x1);
xlabel('$n$', 'interpreter', 'latex', 'FontSize', 13);
ylabel('$h*\tilde w-\tilde x$', 'interpreter', 'latex', 'FontSize', 13);
f = gcf; f.Position(3:4) = 7/10 * [560 420];
ylim(YLimX);


%% diff autocorr

% figure;
% plot(xcorr(x1));
% YLimR = max(abs(get(gca, 'YLim') )) * [-1 1];
% ylim(YLimR);
% 
% figure;
% plot(xcorr(x2));
% ylim(YLimR);
% 
% figure;
% plot(xcorr(x2)-xcorr(x1));
% ylim(YLimR);