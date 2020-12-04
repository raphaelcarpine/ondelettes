n = 10000;

U = ones(1, n);
N = 0.1 * exp(2i*pi*rand(1, n));

R1 = abs(abs(sum((U+N).^2)) - n)
R2 = abs(sum(abs(U+N).^2) - n)