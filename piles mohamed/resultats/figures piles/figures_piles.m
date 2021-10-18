saveFolder = 'piles mohamed\resultats\figures piles\figures';

capts = {'acc1x', 'acc2x', 'acc3x', 'acc4x', 'acc5x'};
nCapts = [4 5 4 5 5 5 4 5 4 5];

for pile = 1:10
    nCapt = nCapts(pile);
    
    % sans sable
    fctDefModale = defModales(pile, capts(1:nCapt), '', nan);
    fig = fctDefModale(nan(1, nCapt));
    fig.Position(4) = 500;
    set(fig,'renderer','Painters');
    saveas(fig, fullfile(saveFolder, sprintf('pile%u', pile)), 'epsc');
    
    % avec sable
    fctDefModaleSable = defModales(pile, capts(1:nCapt), '', 0);
    figSable = fctDefModaleSable(nan(1, nCapt));
    figSable.Position(4) = 500;
    set(figSable,'renderer','Painters');
    saveas(figSable, fullfile(saveFolder, sprintf('pile%u_sable', pile)), 'epsc');
end