% natural freqs 
Fnat = [6, 10];
Znat = 0.01 * [1, 0.8];
Ampl_nat = [1, 0.8; 0.8, 1.3];
Phi_nat = [1, 3; -3, 2];

% train
Ftrain = 4;
Freqs_train = Ftrain * (0:4);
Ampl_train = [0., 0.8, 0.5, 0.4, 0.2];
Phi_train = [0, 3, -1, -1, -2];

% times
t0 = 0;
tf = 10;
t1 = 2;
t2 = 6;
t = linspace(t0, tf, 10000);

% signal
x_train = Ampl_train' .* cos(2*pi* Freqs_train' * t + Phi_train') .* (t>=t1 & t<=t2);
x_nat_1 = Ampl_nat(1, :)' .* cos(2*pi* Fnat' * (t-t1) + Phi_nat(1, :)') .* exp(-2*pi * (Znat'.*Fnat') * (t-t1)) .* (t>=t1);
x_nat_2 = Ampl_nat(2, :)' .* cos(2*pi* Fnat' * (t-t2) + Phi_nat(2, :)') .* exp(-2*pi * (Znat'.*Fnat') * (t-t2)) .* (t>=t2);

x = sum(x_train) + sum(x_nat_1) + sum(x_nat_2);



% figures
figure;
plot(t, x_train);
hold on
plot(t, x_nat_1);
plot(t, x_nat_2);

figure;
plt = plot(t, x + 0.1*randn(size(x)));



% wavelet
fmin = 2;
fmax = 18;
Q = 15;
MaxRidges = 8;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', inf);


