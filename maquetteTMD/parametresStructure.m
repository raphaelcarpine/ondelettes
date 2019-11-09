load('maquetteTMD/donnees3/mData.mat');



%% frequence propre strucutre

mesures = {structure_001, structure_002, structure_003, structure_004, structure_005};
f0 = nan(1, 5);
for k = 1:5
    mesure = transpose(mesures{k});
    t = (0:(length(mesure)-1)) / freq_ech;
    a = mesure;
    f0(k) = getFreq(t, a, [6, 7], 0.001);
end




%% frequence propre strucutre masse ajoutée

mesures2 = {structure_masse_ajoutee_001, structure_masse_ajoutee_002, structure_masse_ajoutee_003,...
    structure_masse_ajoutee_004, structure_masse_ajoutee_005};
f02 = nan(1, 5);
for k = 1:5
    mesure = transpose(mesures2{k});
    t = (0:(length(mesure)-1)) / freq_ech;
    a = mesure;
    f02(k) = getFreq(t, a, [6, 7], 0.001);
end


%% calculs

Df0 = std(f0)/sqrt(length(f0));
Df0 = 0.005/sqrt(length(f0));
f0 = mean(f0);

Df02 = std(f02)/sqrt(length(f02));
Df02 = 0.005/sqrt(length(f02));
f02 = mean(f02);

mu = delta_m/((f0/f02).^2-1);

Dmu = mu * sqrt((Sigma_delta_m/delta_m)^2 + ( (2*Df0*f0/f02^2)^2 +  (2*Df02*f0/f02^2)^2   )/((f0/f02)^2-1)^2);





disp(['f0 = ', num2str(f0), '+-', num2str(Df0)]);
disp(['mu = ', num2str(mu), '+-', num2str(Dmu)]);




























