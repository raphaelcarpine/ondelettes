load('mur silvia\modesFourier\ModesLMS.mat');


load('mur silvia\modesOndelette\allModalQuantities.mat');

allFreqs = AllModalQuantities.freqs;
allShapes = AllModalQuantities.shapes;
allDamps = AllModalQuantities.damps;

%%

P = [0, 6, 7];

tableMat = 0;
lineMat = 1;
columnMat = 1;


%%

for indp = 1:3
    p = P(indp);
    
    for mode = 1:nbModes(indp)
        
        % nb transient
        nbTransients = length(allFreqs{indp}{mode});
        tableMat(lineMat, columnMat) = nbTransients;
        
        
        % freq
        meanFreqs{indp}(mode) = mean(allFreqs{indp}{mode});
        stdFreqs{indp}(mode) = std(allFreqs{indp}{mode});
        
        
        % damping
        meanDamps{indp}(mode) = mean(allDamps{indp}{mode});
        stdDamps{indp}(mode) = std(allDamps{indp}{mode});
        
        
        % shape
        meanShapes{indp}(mode, :) = mean(allShapes{indp}{mode}, 1);
        stdShapes{indp}(mode, :) = std( real( allShapes{indp}{mode}), 0, 1) + 1i*std( imag( allShapes{indp}{mode}), 0, 1);
        
        
        meanFreq = meanFreqs{indp}(mode);
        meanShape = transpose(meanShapes{indp}(mode, :));
        zeta = meanDamps{indp}(mode);
        
        % erreur stat
        errorFreq = stdFreqs{indp}(mode) / sqrt(nbTransients{indp}(mode));
        errorDamp = stdDamps{indp}(mode) / sqrt(nbTransients{indp}(mode));
        errorShape = norm( stdShapes{indp}(mode, :)) / sqrt(nbTransients{indp}(mode));
        errorShape = errorShape / norm(meanShape);
        shapeI = norm( imag( meanShape));
        errorShapeI = norm( imag( stdShapes{indp}(mode, :))) / sqrt(nbTransients{indp}(mode));
        
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
            disp(['freq : ', num2str(meanFreq), ' ; error : ', num2str(100*errorFreq/meanFreq), '%',...
                ' ; freq lms : ', num2str(ModesLMS(indp, mode).freq), ' ; error lms : ', num2str(100*errorFreqF), '%']);
            disp(['MAC : ', num2str(100*mac), '% ; ', ' ; shape error : ', num2str(100*errorShape), '%',...
                ' ; shape error lms : ', num2str(100*errorShapeF), '%']);
            disp(['imaginary shape : ', num2str(shapeI), ' ; ', 'imaginary shape error : ',...
                num2str(100*errorShapeI/shapeI), '%']);
            disp(['amort : ', num2str(100*zeta), '% ; error : ', num2str(100*errorDamp), '%',...
                ' ; relative error : ', num2str(100*errorDamp/zeta), '%',...
                ' ; amort lms : ', num2str(100*ModesLMS(indp, mode).damping), '%']);
            
            disp(['I : ', num2str(100*nonPropIndex(meanShape)), '%',...
                ' ; I lms : ', num2str(100*nonPropIndex(ModesLMS(indp, mode).shape)), '%']);
        end
        
        if plotGlobal
            title = ['P', num2str(p), 'M', num2str(mode), '_freq=', num2str(meanFreq), '_damping=', num2str(100*zeta)];
            fig = plotModShape(real(meanShape), title);
            
            title2 = [title, '_complex'];
            fig2 = plotComplexModShape(meanShape, title2);
            
            % enregistrement
            if saveFigs
                directory = 'mur silvia\modesOndelette\save\';
                savefig(fig, [directory, title, '.fig']);
                saveas(fig, [directory, title, '.png']);
                savefig(fig2, [directory, title2, '.fig']);
                saveas(fig2, [directory, title2, '.png']);
            end
        end
        
    end
end
