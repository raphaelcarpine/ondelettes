%% info

%29277:ch1,29279:ch1,29279:ch2,29279:ch3,29277:ch2,29277:ch3,29278:ch1,29278:ch2,29278:ch3,29280:ch1,29280:ch2,29280:ch3,29281:ch1,29281:ch2,29281:ch3,29282:ch1,29282:ch2,29282:ch3,29283:ch1,29283:ch2,29283:ch3,29284:ch1,29284:ch2,29284:ch3
ordre = [1.1, 3.1, 3.2, 3.3, 1.2, 1.3, 2.1, 2.2, 2.3, 4.1, 4.2, 4.3, 5.1, 5.2, 5.3, 6.1, 6.2, 6.3, 7.1, 7.2, 7.3, 8.1, 8.2, 8.3]; % 2.3 : capteur 2 ch3
[~, ordre] = sort(ordre);


%%

Data = readmatrix('tests capteurs/SensorConnectData3.csv');

Data = Data(:, 2:end);
Data = transpose(Data);

kti = floor(size(Data, 2)/2);
while kti >=1 && ~any(isnan(Data(:, kti)))
    kti = kti-1;
end
kti = kti+1;
ktf = floor(size(Data, 2)/2);
while ktf <= size(Data, 2) && ~any(isnan(Data(:, ktf)))
    ktf = ktf+1;
end
ktf = ktf-1;
Data = Data(:, kti:ktf);
% Data = Data(ordre, :);

% Data = Data(:, 1:floor(end/4));
% Data = Data(3:3:end, :);

Fe = 256;
Fc = 104;
Fband = Fe/2;
dt = 1/Fe;

t = dt * (1:size(Data, 2));

figure;
plt = plot(t, Data);

WaveletMenu('WaveletPlot', plt);

%% découpage pour moyennage

meanData = mean(Data, 2);
Data = Data - meanData;
n = size(Data, 2);
N = 50;
n0 = floor(n/N);
Data0 = cell(1, N);
for k = 1:N
    Data0{k} = Data(:, (k-1)*n0+1:k*n0);
end


%% calcul densité de puissance
PSD = 0;
for k = 1:N
    FTdata = fft(Data0{k}, [], 2);
    FTdata = FTdata(:, 1:floor(n0/2+1));
    PSDk = 1/(n0*Fe) * abs(FTdata).^2;
    PSDk(:, 2:end) = 2*PSDk(:, 2:end);
    PSD = PSD + PSDk;
end
PSD = PSD/N;

figure;
plot(1/(dt*n0) * (0:size(PSD, 2)-1), sqrt(PSD)/9.81*1e6);
xlabel('Fréquence [Hz]');
ylabel('Densité spectrale de bruit [µg/sqrt(Hz)]');

P = mean(Data.^2, 2);
B = sqrt(P / Fband);
for b = B
    fprintf('%.f µg/sqrtHz\n', b/9.81*1e6);
end


%% calcul norme

for capteur = 1:8
    normCapteur = norm(meanData((capteur-1)*3+1:capteur*3));
    fprintf('capteur %u : norme = %.3f m/s² (%.1f%% error)\n', [capteur, normCapteur, 100*(normCapteur/9.81-1)]);
end

