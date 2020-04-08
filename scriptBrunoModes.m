% choix precision
ct = 3;
cf = 5;


% construction du signal
f1 = 10;
f2 = 12;
zeta1 = 0.01;
zeta2 = 0.01;
shape1 = [1; 3; -1; 2];
shape2 = [1; -2; -1; 1+0.1i];
T = 10;

t = linspace(0, T, 10000);
X = shape1 * exp(2i*pi*f1*t - zeta1*2*pi*f1*t);
X = X + shape2 * exp(2i*pi*f2*t - zeta2*2*pi*f2*t);
X = real(X);

% choix Q
Dt1 = 1/(zeta1*2*pi*f1);
Df1 = f2-f1;
[Qmin, Qmax, Qz] = getBoundsQ(f1, Df1, Dt1, T, ct, cf);
disp('premier mode :');
disp(['Qmin = ', num2str(Qmin), ' ; Qmax = ', num2str(Qmax), ' ; Qz = ', num2str(Qz)]);

Dt2 = 1/(zeta2*2*pi*f2);
Df2 = f2-f1;
[Qmin, Qmax, Qz] = getBoundsQ(f2, Df2, Dt2, T, ct, cf);
disp('deuxième mode :');
disp(['Qmin = ', num2str(Qmin), ' ; Qmax = ', num2str(Qmax), ' ; Qz = ', num2str(Qz)]);

Q = 30;
disp(' ');
disp(['Q = ', num2str(Q)]);
disp(' ');

% CWT
NbMaxRidges = 2;
NbMaxParallelRidges = inf;
fmin = 8;
fmax = 114;
nbFreqs = 300;

[time, freq, shapes, amplitudes] = getModesSingleRidge(t, X, Q, fmin, fmax, nbFreqs,... % attention t et time sont différents ! time est le vecteur de temps du ridge uniquement
    'NbMaxRidges', NbMaxRidges, 'NbMaxParallelRidges', NbMaxParallelRidges);


% affichage
for mode = 1:length(time)
    % partie imaginaire deformee
    figure;
    plot(time{mode}, imag(shapes{mode}));
    xlabel('time');
    ylabel('Im\phi');
    
    % partie reelle deformee
    figure;
    plot(time{mode}, real(shapes{mode}));
    xlabel('time');
    ylabel('Re\phi');
    
    % amplitude
    figure;
    plot(time{mode}, abs(amplitudes{mode}));
    set(gca, 'YScale', 'log')
    xlabel('time');
    ylabel('ampl');
    
    % frequence
    figure;
    plot(time{mode}, freq{mode});
    xlabel('time');
    ylabel('freq');
    
    % moyennes
    fmoy = mean(freq{mode});
    shapemoy = mean(shapes{mode}, 2);
    disp(['mode ', num2str(mode)]);
    disp(['freq moy : ', num2str(fmoy)]);
    disp('shape moy : ');
    disp(shapemoy);
end
