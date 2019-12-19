directory = 'mur silvia\modesFourier\';

files{1} = 'LMS_ModalShapes_P0.xlsx';
files{2} = 'LMS_ModalShapes_P6.xlsx';
files{3} = 'LMS_ModalShapes_P7.xlsx';

P = [0 6 7];

ModesLMS = struct([]);

for k = 1:3
    tab = xlsread([directory, files{k}]);
    
    p = P(k);
    
    for mode = 1:3
        L0 = 11*(mode-1);
        freq = tab(L0 + 1, 2);
        
        shape = tab((L0+2):(L0+10), 3) + 1i*tab((L0+2):(L0+10), 4);
        
        shape = shape / sqrt(transpose(shape) * shape);
        shape = sign(real(shape(1))) * shape;
        
        % enregistrement
        ModesLMS(k, mode).freq = freq;
        ModesLMS(k, mode).shape = shape;
    end
end


save([directory, 'ModesLMS.mat'], 'ModesLMS');

