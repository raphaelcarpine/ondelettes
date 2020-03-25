load('mur silvia\modesFourier\ModesLMS.mat');
load('mur silvia\modesOndelette\allModalQuantities.mat');
load('mur silvia\modesOndelette2\ModesWvlt2.mat');

allShapes = AllModalQuantities.shapes;

%%
P = [0, 6, 7];
nbModes = [3, 4, 3];


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
    disp('LMS CWT')
    disp(MAC12);
    disp('LMS CWT2')
    disp(MAC13);
    disp('CWT CWT2')
    disp(MAC23);
    disp([newline, newline']);
    
    cmap = gray(256);
    cmap = [cmap; cmap(256:-1:1, :)];
    cmap(1, :) = [1, 0, 0];
    
    MAC12(isnan(MAC12)) = -1;
    MAC13(isnan(MAC13)) = -1;
    MAC23(isnan(MAC23)) = -1;
    
    fig = figure;
    ax = axes(fig);
    imagesc(ax, MAC12, [-1, 1]);
    colormap(ax, cmap);
    ylabel(ax, 'LMS');
    xlabel(ax, 'CWT');
    title(ax, ['P', num2str(p)]);
    
    fig = figure;
    ax = axes(fig);
    imagesc(ax, MAC13, [-1, 1]);
    colormap(ax, cmap);
    ylabel(ax, 'LMS');
    xlabel(ax, 'CWT2');
    title(ax, ['P', num2str(p)]);
    
    fig = figure;
    ax = axes(fig);
    imagesc(ax, MAC23, [-1, 1]);
    colormap(ax, cmap);
    ylabel(ax, 'CWT');
    xlabel(ax, 'CWT2');
    title(ax, ['P', num2str(p)]);
end




