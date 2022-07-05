[T, X, dt, chDist] = getDataRadio('Mesure 5', 'mesures preliminaires');

selectDists = [];%[21.8, 31.6, 45.9, 56.4];%[31.6];
if ~isempty(selectDists)
    Kdists = ismember(chDist, selectDists);
    chDist = chDist(Kdists);
    X = X(Kdists, :);
end

figure;
plt = plot(T, X);
xlabel('Temps [s]');
ylabel('Déplacement vertical [mm]');
legend(strsplit(num2str(chDist)));
selectLine();


WaveletMenu('WaveletPlot', plt);