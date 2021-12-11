function modesTable = groupShapesByMac(Shapes)
%GROUPSHAPESBYMAC Summary of this function goes here
%   Shapes{k_pont}(:, k_freq)

limitMAC = 0.5; % en dessous de ce MAC, modes différents
realMAC = true; % prendre le mac des parties réelles

N = size(Shapes{1}, 1); % nb de ddl

modesTableCell = {}; % groupes de fréquences
for kf = 1:size(Shapes{1}, 2)
    modesTableCell{end+1} = [1; kf];
end

for kp = 2:length(Shapes)
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
    [I, MACs] = matchShapes(modesTableMeanShapes, Shapes{kp}, realMAC);
    for kf = 1:size(Shapes{kp}, 2)
        if isnan(MACs(kf)) || MACs(kf) < limitMAC
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





function [I, MACs] = matchShapes(shapes1, shapes2, realMAC)

if size(shapes1, 2) >= size(shapes2, 2)
    shapesL = shapes1;
    shapesS = shapes2;
else
    shapesL = shapes2;
    shapesS = shapes1;
end

if realMAC
    shapesL = real(shapesL);
    shapesS = real(shapesS);
end

MACmat = (shapesL'*shapesS).^2 ./ (diag(shapesL'*shapesL) * diag(shapesS'*shapesS).');

[~, I] = max(MACmat);

if ~ all(diff(sort(I)) > 0) % tous les éléments de I ne sont pas distincts
    % séparation dans l'ordre
    for kf = 2:length(I)
        MACmatk = MACmat(:, kf);
        MACmatk(I(1:kf-1)) = 0;
        [~, I(kf)] = max(MACmatk);
    end
    
    % test réarrangement
    chgt = true;
    while chgt
        chgt = false;
        for kch1 = 1:length(I)
            for kch2 = kch1+1:length(I)
                sumMAC0 = sum(diag(MACmat(I, 1:length(I))));
                I2 = I;
                I2([kch1, kch2]) = I2([kch2, kch1]);
                sumMACch = sum(diag(MACmat(I2, 1:length(I2))));
                if sumMACch > sumMAC0
                    I = I2;
                    chgt = true;
                    break
                end
            end
            if chgt
                break
            end
        end
    end
end

% calcul macs
MACs = diag(MACmat(I, 1:length(I)));


if size(shapes1, 2) < size(shapes2, 2) % inversion pour shape2 dans shape1
    I0 = I;
    MACs0 = MACs;
    I = nan(1, size(shapes2, 2));
    MACs = nan(1, size(shapes2, 2));
    for ki0 = 1:length(I0)
        I(I0(ki0)) = ki0;
        MACs(I0(ki0)) = MACs0(ki0);
    end
end

end

