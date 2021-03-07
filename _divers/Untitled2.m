k1 = 1;
k2 = 0.5;
m = 1;
x0 = 4;

g = @(x) k1*x + (k2-k1)*(x-x0).*(x>=x0) + (k2-k1)*(x+x0).*(x<=-x0);
G = @(x) k1*x.^2/2 + (k2-k1)*(x-x0).^2.*(x>=x0)/2 + (k2-k1)*(x+x0).^2.*(x<=-x0)/2;
I = @(A) sqrt(m/k1)*asin(x0./sqrt(A.^2+(k2-k1)/k1*(A-x0).^2)) + sqrt(m/k2)*acos(k1*x0./(k2*A+(k1-k2)*x0));


x = linspace(-10, 10, 1000);
A = linspace(0, 10, 1000);

%%
figure;
plot(x, g(x));

%%
G2 = nan(size(x));
for ix = 1:length(x)
    G2(ix) = integral(g, 0, x(ix));
end

figure;
plot(x, G(x));
hold on
plot(x, G2, '--');

%%
I2 = nan(size(A));
for ia = 1:length(A)
    I2(ia) = integral(@(x) 1./sqrt(2*(G(A(ia))-G(x))/m), 0, A(ia));
end

figure;
plot(A, I(A));
hold on
plot(A, I2, '--');