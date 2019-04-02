omega0 = '2*pi';
epsilon = 0.1;
n = 1;
delta = 0.05;

d2x = '-omega0^2*x-(abs(x)<=delta)*epsilon*sign(dx)*abs(dx)^n';

T = 100;


systemeQuelconque({'x'}, {d2x}, {'omega0', 'epsilon', 'n', 'delta'}, {omega0, epsilon, n, delta}, 0, 1, 'T', T);

RegressionMenu;