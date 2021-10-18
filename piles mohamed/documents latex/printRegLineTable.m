dataFolder = 'piles mohamed\resultats\reg lin';

pile = 10;
direction = 'x';
g = 8;
mode = 2;

%%

fileName = sprintf('pile%u_%c_%ug_mode%u_reglin.xlsx', [pile, int8(direction), g, mode]);
A = readmatrix(fullfile(dataFolder, fileName));

A(:, 4:5) = 100*A(:, 4:5); % pourcentage

for kl = 1:size(A, 1)
    if ~isnan(A(kl, 2))
        fprintf(' & %i & & %.1f & %.1f & & %.1f & %.1f \\\\ \n', A(kl, :));
    else
        fprintf(' & %i & & / & / & & / & / \\\\ \n', A(kl, 1));
    end
end