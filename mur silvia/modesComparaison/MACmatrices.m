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
    
    for modeX = 1:nbModes(indp)
        for modeY = 1:nbModes(indp)
            
            shapeF = ModesLMS(indp,modeX).shape;
            meanShapeX = transpose( mean( allShapes{indp}{modeX}, 1));
            meanShapeY = transpose( mean( allShapes{indp}{modeY}, 1));
            if modeY <= size(ModesWvlt2, 2) && ~isempty(ModesWvlt2(indp,modeY).shape)
                shapeCWT2 = ModesWvlt2(indp,modeY).shape;
            else
                shapeCWT2 = nan(9, 1);
            end
            
            % mac 12
            mac = abs(meanShapeY'*shapeF)^2 / ((meanShapeY'*meanShapeY) * (shapeF'*shapeF));
            MAC12(modeX, modeY) = mac;
            
            % mac 13
            mac = abs(shapeCWT2'*shapeF)^2 / ((shapeCWT2'*shapeCWT2) * (shapeF'*shapeF));
            MAC13(modeX, modeY) = mac;
            
            % mac 23
            mac = abs(shapeCWT2'*meanShapeX)^2 / ((shapeCWT2'*shapeCWT2) * (meanShapeX'*meanShapeX));
            MAC23(modeX, modeY) = mac;
        end
    end
    
    disp(['%%%%%% P', num2str(p), ' %%%%%%', newline]);
    disp('LSCF CWT')
    disp(100*MAC12);
    disp('LSCF CWT2')
    disp(100*MAC13);
    disp('CWT CWT2')
    disp(100*MAC23);
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


