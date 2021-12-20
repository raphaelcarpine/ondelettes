function modesTable = groupModShapes0(Freqs, Shapes)
%GROUPSHAPESBYMAC Summary of this function goes here
%   Shapes{k_pont}(:, k_freq)
%   Freqs{k_pont}(k_freq)

limitFreq = 0.1; % difference relative, au dessus, modes différents
limitMAC = 0.; % en dessous de ce MAC, modes différents
realMAC = true; % prendre le mac des parties réelles

N = size(Shapes{1}, 1); % nb de ddl

modesTableCell = {}; % groupes de fréquences
for kf = 1:size(Shapes{1}, 2)
    modesTableCell{end+1} = [1; kf];
end

for kp = 2:length(Shapes)
    % construction des freqs moyennes
    modesTableMeanFreqs = zeros(1, length(modesTableCell));
    for kmod = 1:length(modesTableCell)
        for kp2 = 1:size(modesTableCell{kmod}, 2)
            modesTableMeanFreqs(kmod) = modesTableMeanFreqs(kmod) +...
                Freqs{modesTableCell{kmod}(1, kp2)}(modesTableCell{kmod}(2, kp2));
        end
        modesTableMeanFreqs(kmod) = modesTableMeanFreqs(kmod) / size(modesTableCell{kmod}, 2);
    end
    
    % construction deformees moyennes
    modesTableMeanShapes = zeros(N, length(modesTableCell));
    for kmod = 1:length(modesTableCell)
        for kp2 = 1:size(modesTableCell{kmod}, 2)
            modesTableMeanShapes(:, kmod) = modesTableMeanShapes(:, kmod) +...
                Shapes{modesTableCell{kmod}(1, kp2)}(:, modesTableCell{kmod}(2, kp2));
        end
        modesTableMeanShapes(:, kmod) = modesTableMeanShapes(:, kmod) / size(modesTableCell{kmod}, 2);
    end
    
    % match des deformees avec les moyennes
    [I, Dfreqs, MACs] = matchModes(modesTableMeanFreqs, Freqs{kp},...
        modesTableMeanShapes, Shapes{kp}, realMAC);
    for kf = 1:size(Shapes{kp}, 2)
        if isnan(Dfreqs(kf)) || Dfreqs(kf) > limitFreq || MACs(kf) < limitMAC
            modesTableCell{end+1} = [kp; kf];
        else
            modesTableCell{I(kf)}(:, end+1) = [kp; kf];
        end
    end
end

% construction tableau modes
modesTable = nan(length(modesTableCell), length(Shapes));
for km = 1:length(modesTableCell)
    for ks = 1:size(modesTableCell{km}, 2)
        modesTable(km, modesTableCell{km}(1, ks)) = modesTableCell{km}(2, ks);
    end
end

end





function [I, Dfreqs, MACs] = matchModes(freqs1, freqs2, shapes1, shapes2, limitFreq, limitMAC, realMAC)

if length(freqs1) >= length(freqs2)
    freqsL = freqs1;
    freqsS = freqs2;
    shapesL = shapes1;
    shapesS = shapes2;
else
    freqsL = freqs2;
    freqsS = freqs1;
    shapesL = shapes2;
    shapesS = shapes1;
end

if realMAC
    shapesL = real(shapesL);
    shapesS = real(shapesS);
end

MACmat = (shapesL'*shapesS).^2 ./ (diag(shapesL'*shapesL) * diag(shapesS'*shapesS).');
Dfreqsmat = 2 * abs(freqsL.' - freqsS) ./ (freqsL.' + freqsS);

% [~, I] = max(MACmat);
[~, I]= min(Dfreqsmat);

if ~ all(diff(sort(I)) > 0) % tous les éléments de I ne sont pas distincts
    % séparation dans l'ordre
    for kf = 2:length(I)
        Dfreqsmatk = Dfreqsmat(:, kf);
        Dfreqsmatk(I(1:kf-1)) = 2;
        [~, I(kf)] = min(Dfreqsmatk);
    end
end

% test réarrangement MAC
chgt = true;
while chgt
    chgt = false;
    for kch1 = 2:length(I)-1
        % echange avant
        kch21 = kch1 - 1;
        sumMAC0 = sum(diag(MACmat(I, 1:length(I))));
        I21 = I;
        I21([kch1, kch21]) = I21([kch21, kch1]);
        sumMACch1 = sum(diag(MACmat(I21, 1:length(I21))));
        %echange après
        kch22 = kch1 + 1;
        I22 = I;
        I22([kch1, kch22]) = I22([kch22, kch1]);
        sumMACch2 = sum(diag(MACmat(I22, 1:length(I22))));
        if sumMACch1 >= sumMACch2 && sumMACch1 > sumMAC0
            I = I21;
            chgt = true;
            break
        elseif sumMACch2 > sumMAC0
            I = I22;
            chgt = true;
            break
        end
        if chgt
            break
        end
    end
end

% calcul macs
MACs = diag(MACmat(I, 1:length(I)));
Dfreqs = diag(Dfreqsmat(I, 1:length(I)));


if size(shapes1, 2) < size(shapes2, 2) % inversion pour shape2 dans shape1
    I0 = I;
    Dfreqs0 = Dfreqs;
    MACs0 = MACs;
    I = nan(1, length(freqs2));
    Dfreqs = nan(1, length(freqs2));
    MACs = nan(1, length(freqs2));
    for ki0 = 1:length(I0)
        I(I0(ki0)) = ki0;
        Dfreqs(I0(ki0)) = Dfreqs0(ki0);
        MACs(I0(ki0)) = MACs0(ki0);
    end
end

end

