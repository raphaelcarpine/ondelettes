m = 1;
k = (2*pi*10)^2;
c = 2*0.01*2*pi*10;

dt = 0.001;
T = 20;
t = 0:dt:T;
f = 1*randn(size(t));

D = @(t, X) [X(2);
    1/m * (f(floor(t/dt)+1) - k*X(1) - c*X(2))];

X0 = [0; 1];
X = RK4(D, t, X0);
x = X(1, :);

figure;
plot(t, x);
xlabel('Temps');
ylabel('Déplacement');
xticks([0 T]);
xticklabels({'0', 'T'});
yticks([]);