f0 = 2;

syst = systLin(1, (2*pi*f0)^2, 0.01*2*pi*f0);

dt = 0.01;
T = 3600;
t = 0:dt:T;
f = randn(size(t));
n = 1;
while 1
    n = n + round(exprnd(30)/dt);
    if n > length(f)
        break
    end
    f(n) = f(n) + 100;
end
x = syst.response(f, dt, 1);

% flag debut choc
Dt = 0.5;
Dt = ones(1, round(Dt/dt));
Dt = Dt / length(Dt);

flag = conv(f, Dt);
flag = flag(1:length(f));

%%
figure;
yyaxis left
plot(t, f);
plot(t, flag);
yyaxis right
plt = plot(t, x);

%%

MotherWavelet = 'morlet';
Q = 6;
fmin = 1;
fmax = 3;

ridge = RidgeExtract(t, x, Q, fmin, fmax, 100, 'MotherWavelet', MotherWavelet,...
        'NbMaxRidges', 1, 'MaxSlopeRidge', inf, 'NbMaxParallelRidges', inf);

    
    