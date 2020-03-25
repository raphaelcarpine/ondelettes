load('mur silvia\modesFourier\ModesLMS.mat');

load('mur silvia\modesOndelette\allModalQuantities.mat');

allFreqs = AllModalQuantities.freqs;
allShapes = AllModalQuantities.shapes;
allDamps = AllModalQuantities.damps;

%%
verb = true;
plotGlobal = true;
saveFigs = true;

directory = 'mur silvia\modesOndelette\save\';

%%

meanFreqs = {[], [], []};
meanShapes = {[], [], []};
meanDamps = {[], [], []};
nbTransients = {[], [], []};
stdFreqs = {[], [], []};
stdShapes = {[], [], []};
stdDamps = {[], [], []};

P = [0, 6, 7];


%%

if saveFigs
    delete([directory, '*']);
end

for indp = 1:3
    p = P(indp);
    
    if verb
        disp(' ');
        disp(' ');
        disp(['~~~~~ P', num2str(p), ' : mean, std']);
    end
    
    for mode = 1:nbModes(indp)
        meanFreqs{indp}(mode) = mean(allFreqs{indp}{mode});
        stdFreqs{indp}(mode) = std(allFreqs{indp}{mode});
        meanDamps{indp}(mode) = mean(allDamps{indp}{mode});
        stdDamps{indp}(mode) = std(allDamps{indp}{mode});
        
        meanShapes{indp}(mode, :) = mean(allShapes{indp}{mode}, 1);
        stdShapes{indp}(mode, :) = std( real( allShapes{indp}{mode}), 0, 1) + 1i*std( imag( allShapes{indp}{mode}), 0, 1);
        
        nbTransients{indp}(mode) = length(allFreqs{indp}{mode});
        
        meanFreq = meanFreqs{indp}(mode);
        meanShape = transpose(meanShapes{indp}(mode, :));
        zeta = meanDamps{indp}(mode);
        
        % erreur stat
        stdFreq = stdFreqs{indp}(mode);
        stdDamp = stdDamps{indp}(mode);
        stdShape = norm( stdShapes{indp}(mode, :));
        stdShape = stdShape / norm(meanShape);
        shapeI = norm( imag( meanShape));
        stdShapeI = norm( imag( stdShapes{indp}(mode, :)));
        
        % comparaison fourier
        errorFreqF = abs(meanFreq - ModesLMS(indp, mode).freq) / ModesLMS(indp, mode).freq;
        
        shapeF = ModesLMS(indp, mode).shape;
        errorShapeF = meanShape - shapeF;
        errorShapeF = sqrt( errorShapeF'*errorShapeF / (shapeF'*shapeF));
        
        % modal assurance criterion
        mac = abs(meanShape'*shapeF)^2 / ((meanShape'*meanShape) * (shapeF'*shapeF));
        
        if verb
            disp(' ');
            disp(['~ mode', num2str(mode), ' (', num2str(nbTransients{indp}(mode)), ' transients)']);
            disp(['freq : ', num2str(meanFreq), ' ; std : ', num2str(stdFreq),...
                ' ; freq lms : ', num2str(ModesLMS(indp, mode).freq), ' ; error lms : ', num2str(100*errorFreqF), '%']);
            disp(['MAC : ', num2str(100*mac), '% ; ', ' ; shape cov : ', num2str(100*stdShape), '%',...
                ' ; shape error lms : ', num2str(100*errorShapeF), '%']);
            disp(['imaginary shape : ', num2str(shapeI), ' ; ', 'imaginary shape cov : ',...
                num2str(100*stdShapeI/shapeI), '%']);
            disp(['amort : ', num2str(100*zeta), '% ; std : ', num2str(100*stdDamp), '%',...
                ' ; cov : ', num2str(100*stdDamp/zeta), '%',...
                ' ; amort lms : ', num2str(100*ModesLMS(indp, mode).damping), '%']);
            
            disp(['I : ', num2str(100*nonPropIndex(meanShape)), '%',...
                ' ; I lms : ', num2str(100*nonPropIndex(ModesLMS(indp, mode).shape)), '%']);
        end
        
        if plotGlobal
            title = ['P', num2str(p), 'M', num2str(mode), 'CWT'];
            fig = plotModShape(real(meanShape), title);
            
            title2 = [title, '_complex'];
            fig2 = plotComplexModShape(meanShape, title2);
            
            % enregistrement
            if saveFigs
                savefig(fig, [directory, title, '.fig']);
                saveas(fig, [directory, title, '.eps'], 'epsc');
                saveas(fig, [directory, title, '.png']);
                savefig(fig2, [directory, title2, '.fig']);
                saveas(fig2, [directory, title2, '.eps'], 'epsc');
                saveas(fig2, [directory, title2, '.png']);
            end
        end
        
    end
end
