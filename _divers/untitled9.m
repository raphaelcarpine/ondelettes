n = 100000;
x = rand(1, n);
y = 3-0.5*x + randn(size(x));
t = linspace(0, 10, length(x));

[X, Y, stdY, K, Xlims] = averagingScatter3(x, y, 0.1, 'log', t, 1);

%%
figure;
scatter(x, y);

figure;
errorbar(X, Y, 1.96*stdY./sqrt(K-1));

figure;
bar(X, K);