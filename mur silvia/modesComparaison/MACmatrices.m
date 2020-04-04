load('mur silvia\modesFourier\ModesLMS.mat');
load('mur silvia\modesOndelette\allModalQuantities.mat');
load('mur silvia\modesOndelette2\ModesWvlt2.mat');

allShapes = AllModalQuantities.shapes;

%%
P = [0, 6, 7];
nbModes = [3, 4, 3];

saveFigs = true;
directory = 'mur silvia\modesComparaison\save\';


for indp = 1:3
    p = P(indp);
    
    MAC12 = nan(nbModes(indp));
    MAC13 = nan(nbModes(indp));
    MAC23 = nan(nbModes(indp));
    
    MAC11 = nan(nbModes(indp));
    MAC22 = nan(nbModes(indp));
    MAC33 = nan(nbModes(indp));
    
    for modeX = 1:nbModes(indp)
        for modeY = 1:nbModes(indp)
            
            shapeLSCFx = ModesLMS(indp, modeX).shape;
            shapeLSCFy = ModesLMS(indp, modeY).shape;
            shapeCWT1x = transpose( mean( allShapes{indp}{modeX}, 1));
            shapeCWT1y = transpose( mean( allShapes{indp}{modeY}, 1));
            if modeX <= size(ModesWvlt2, 2) && ~isempty(ModesWvlt2(indp, modeX).shape)
                shapeCWT2x = ModesWvlt2(indp, modeX).shape;
            else
                shapeCWT2x = nan(9, 1);
            end
            if modeY <= size(ModesWvlt2, 2) && ~isempty(ModesWvlt2(indp, modeY).shape)
                shapeCWT2y = ModesWvlt2(indp, modeY).shape;
            else
                shapeCWT2y = nan(9, 1);
            end
            
            % MACs
            MAC12(modeX, modeY) = MACnb(shapeLSCFx, shapeCWT1y);
            MAC13(modeX, modeY) = MACnb(shapeLSCFx, shapeCWT2y);
            MAC23(modeX, modeY) = MACnb(shapeCWT1x, shapeCWT2y);
            MAC11(modeX, modeY) = MACnb(shapeLSCFx, shapeLSCFy);
            MAC22(modeX, modeY) = MACnb(shapeCWT1x, shapeCWT1y);
            MAC33(modeX, modeY) = MACnb(shapeCWT2x, shapeCWT2y);
            
        end
    end
    
    disp(['%%%%%% P', num2str(p), ' %%%%%%', newline]);
    disp('LSCF CWT')
    disp(100*MAC12);
    disp('LSCF CWT2')
    disp(100*MAC13);
    disp('CWT CWT2')
    disp(100*MAC23);
    disp('LSCF')
    disp(100*MAC11);
    disp('CWT1')
    disp(100*MAC22);
    disp('CWT2')
    disp(100*MAC33);
    disp([newline, newline']);
    
    
    [fig1, ax] = dispMACmat(MAC12);
    ylabel(ax, 'LSCF');
    xlabel(ax, 'CWT');
    title(ax, ['P', num2str(p)]);
    title1 = ['P', num2str(p),'_MAC_LSCF_CWT'];
    set(fig1, 'Name', title1);
    
    [fig2, ax] = dispMACmat(MAC13);
    ylabel(ax, 'LSCF');
    xlabel(ax, 'CWT 2');
    title(ax, ['P', num2str(p)]);
    title2 = ['P', num2str(p),'_MAC_LSCF_CWT2'];
    set(fig2, 'Name', title2);
    
    [fig3, ax] = dispMACmat(MAC23);
    ylabel(ax, 'CWT');
    xlabel(ax, 'CWT 2');
    title(ax, ['P', num2str(p)]);
    title3 = ['P', num2str(p),'_MAC_CWT_CWT2'];
    set(fig3, 'Name', title3);
    
    if saveFigs
        saveas(fig1, [directory, title1, '.eps'], 'epsc');
        saveas(fig2, [directory, title2, '.eps'], 'epsc');
        saveas(fig3, [directory, title3, '.eps'], 'epsc');
    end
end



function [fig, ax] = dispMACmat(MACmat)

% echelle gris
cmap = gray(256);
cmap = [cmap; cmap(256:-1:1, :)];
cmap(1, :) = [1, 0, 0]; % nan values

% nan
MACmat(isnan(MACmat)) = -1;

% plot
fig = figure;
pos = get(fig, 'Position');
set(fig, 'Position', pos([1 2 4 4]));
ax = axes(fig);
hold(ax, 'on');
imagesc(ax, MACmat, [-1, 1]);
for ix = 1:size(MACmat, 1)
    for iy = 1:size(MACmat, 2)
        if MACmat(ix, iy) ~= -1
            if MACmat(ix, iy) <= 0.7
                textColor = 'black';
            else
                textColor = 'white';
            end
            text(ax, iy, ix, sprintf('%.2f %%', 100*MACmat(ix, iy)), 'Color', textColor,...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end
hold(ax, 'off');
colormap(ax, cmap);
xlim(ax, [0.5, size(MACmat, 1) + 0.5]);
ylim(ax, [0.5, size(MACmat, 2) + 0.5]);
pbaspect(ax, [1 1 1]);
set(ax, 'YDir','reverse')

xticks(ax, 1:size(MACmat, 1));
yticks(ax, 1:size(MACmat, 2));
set(ax, 'TickLength', [0 0])

end


