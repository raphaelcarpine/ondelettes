directory = 'mur silvia\modesFourier\';

files{1} = 'B1_1_P0_03032020.xlsx';
files{2} = 'B1_1_P6_03032020.xlsx';
files{3} = 'B1_1_P7_03032020.xlsx';
nbModes = [3, 4, 3];

% files{1} = 'LMS_ModalShapes_P0.xlsx';
% files{2} = 'LMS_ModalShapes_P6.xlsx';
% files{3} = 'LMS_ModalShapes_P7.xlsx';
% nbModes = [3, 3, 3];

P = [0 6 7];

ModesLMS = struct([]);

for k = 1:3
    tab = xlsread([directory, files{k}]);
    
    p = P(k);
    
    for mode = 1:nbModes(k)
        L0 = 11*(mode-1);
        freq = tab(L0 + 1, 2);
        
        shape = tab((L0+2):(L0+10), 3) + 1i*tab((L0+2):(L0+10), 4);
        
        shape = shape / sqrt(transpose(shape) * shape);
        shape = sign(real(shape(1))) * shape;
        
        try
            damping = tab(L0+2, 5);
        catch
            damping = nan;
        end
        
        % enregistrement
        ModesLMS(k, mode).freq = freq;
        ModesLMS(k, mode).shape = shape;
        ModesLMS(k, mode).damping = damping;
    end
end


save([directory, 'ModesLMS.mat'], 'ModesLMS');

