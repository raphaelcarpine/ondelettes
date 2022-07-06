clear all

%% parametres

testWaveletMenu = 0;
testAffichage = 0;

kradio = 2;

ridgeContinuity = 'none';

projection = 2; % 1: mode 1, 2, 5

% CWT
MotherWavelet = 'morlet';
ct = 3;
ridgeContinuity = 'none'; % 'none', 'simple', 'reverse, 'double', 'slope3'
if length(ridgeContinuity) >= 5 && strcmp(ridgeContinuity(1:5), 'slope')
    slopeTimeConst = str2double(ridgeContinuity(6:end));
    ridgeContinuity = 'slope';
else
    slopeTimeConst = nan;
end
switch projection
    case 1
        fmin = 1;
        fmax = 3;
        Q = 2.5;
    case 2
        fmin = 1;
        fmax = 3;
        Q = 2.5;
    case 5
        fmin = 5.32 - 3;
        fmax = 5.32 + 3;
        Q = 2.5;
    otherwise
        warning('pas de mode sélectionné');
end

%% data

[filePath, fileName, filePathRadio, fileNameRadio, Nlag] = choixData(0);


% acceleros
load(filePath);
A = X.';
T = T.';
[A, T] = removeRedundantData(A, T);
[A, T] = removeNanSignal(A, T);

A = A - mean(A, 2);

switch projection
    case 1
        k1 = find(contains(channelNames, '29280:ch3'));
        k2 = find(contains(channelNames, '40199:ch3'));
        A = A(k1, :) + A(k2, :);
    case 2
        k1 = find(contains(channelNames, '29279:ch3'));
        k2 = find(contains(channelNames, '29281:ch3'));
        k3 = find(contains(channelNames, '40196:ch3'));
        k4 = find(contains(channelNames, '40200:ch3'));
        A = A(k1, :) - A(k2, :) + A(k3, :) - A(k4, :);
    case 5
        k1 = find(contains(channelNames, '29279:ch3'));
        k2 = find(contains(channelNames, '29281:ch3'));
        k3 = find(contains(channelNames, '40196:ch3'));
        k4 = find(contains(channelNames, '40200:ch3'));
        A = A(k1, :) - A(k2, :) - A(k3, :) + A(k4, :);
    otherwise
%         error(' ');
end


% radio
load(filePathRadio);

% restriction intervalle tps radio
ki = 1 + Nlag;
kf = Nlag + length(Tradio);
T = Tradio;
A = A(:, ki:kf);

%% test

if testWaveletMenu
    figure;
    plt = plot(T, A);
    WaveletMenu('WaveletPlot', plt);
    return
end

%% test

if testAffichage
    figure;
    yyaxis left
    plot(T, Xradio(kradio, :));
    xlabel('Temps [s]');
    ylabel('Flèche [mm]');
    yyaxis right
    plot(T, A);
    ylabel('Accélération [m/s²]');
    return
end

%% CWT

% CWT
freqs = linspace(fmin, fmax, 300);
CWT = WvltComp(T, A, freqs, Q, 'MotherWavelet', MotherWavelet, 'DisplayWaitBar', true);


% ridge
ridge = SingleRidgeExtract(T, freqs, CWT, MotherWavelet, Q, ct, ridgeContinuity, slopeTimeConst);
Fridge = ridge.freq;
Aridge = ridge.val;

%% plot

% % plot temporel
% figure;
% plot(ridge.time, ridge.freq);
% xlabel('Temps [s]');
% ylabel('Fréquence [Hz]');

% % plot amplitude
% figure;
% plot(abs(ridge.val), ridge.freq);
% xlabel('Amplitude [m/s²]');
% ylabel('Fréquence [Hz]');

% plot fleche
figure;
plot(Xradio(kradio, T >= ridge.time(1) & T <= ridge.time(end)), ridge.freq);
xlabel('Flèche [mm]');
ylabel('Fréquence [Hz]');


Xtot = Xradio(kradio, T >= ridge.time(1) & T <= ridge.time(end)).';
Ftot = ridge.freq.';
noNaN = ~isnan(Xtot) & ~isnan(Ftot);
Xtot = Xtot(noNaN);
Ftot = Ftot(noNaN);
C = [ones(size(Xtot)), Xtot] \ Ftot;
fprintf('\ncorr = %.0f mHz/mm\n', 1000*C(2));

hold on
plot(get(gca, 'XLim'), C(1) + C(2)*get(gca, 'XLim'), '--r');

% RegressionMenu

%% estimation erreur coeff corr

Ninterv = 10; % nb d'intervalles de division du signal

Nsub = floor(length(Xtot)/Ninterv);
coeffs = nan(1, Ninterv);
for ki = 1:Ninterv
    Xsub = Xtot(1+(ki-1)*Nsub:ki*Nsub);
    Fsub = Ftot(1+(ki-1)*Nsub:ki*Nsub);
    C = [ones(size(Xsub)), Xsub] \ Fsub; % Fsub = [ones(size(Xsub)), Xsub] * C
    coeffs(ki) = C(2);
end

fprintf('\ncorr = %.0f +- %.0f mHz/mm\n', 1000*[mean(coeffs), std(coeffs)/sqrt(Ninterv)]);







