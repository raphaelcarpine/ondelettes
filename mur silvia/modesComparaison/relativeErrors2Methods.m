load('mur silvia\modesFourier\ModesLMS.mat');
load('mur silvia\modesOndelette\allModalQuantities.mat');
load('mur silvia\modesOndelette2\ModesWvlt2.mat');

allFreqs = AllModalQuantities.freqs;
allShapes = AllModalQuantities.shapes;
allDamps = AllModalQuantities.damps;

%%
P = [0, 6, 7];

nbModes = [3, 4, 3];

errorsFreqs = [];
errorsDamps = [];

%%

for indp = 1:3
    p = P(indp);
    
    for mode = 1:nbModes(indp)        
        %% freqs
        modeFreqs = [];
        
        % freq LMS
        modeFreqs = [modeFreqs, ModesLMS(indp,mode).freq];
        
        % mean freqs cwt
        modeFreqs = [modeFreqs, mean(allFreqs{indp}{mode})];
        
        errorsFreqs = [errorsFreqs, diff(modeFreqs)/mean(modeFreqs)];
        
        %% damping
        modeDamps = [];
        
        % damping LMS
        modeDamps = [modeDamps, 100*ModesLMS(indp,mode).damping];
        
        % mean dampings cwt
        modeDamps = [modeDamps, 100*mean(allDamps{indp}{mode})];
        
        errorsDamps = [errorsDamps, diff(modeDamps)/mean(modeDamps)];
        
    end
end

errorsFreqs = abs(errorsFreqs);
errorsDamps = abs(errorsDamps);


disp('~~~~~~~~ Freqs ~~~~~~~~')
disp(' ')
disp(100*errorsFreqs);
disp(['max : ', num2str(max(100*errorsFreqs)), '%']);
disp(['mean : ', num2str(mean(100*errorsFreqs)), '%']);

disp(' ')
disp(' ')
disp('~~~~~~~~ Damps ~~~~~~~~')
disp(' ')
disp(100*errorsDamps);
disp(['max : ', num2str(max(100*errorsDamps)), '%']);
disp(['mean : ', num2str(mean(100*errorsDamps)), '%']);



