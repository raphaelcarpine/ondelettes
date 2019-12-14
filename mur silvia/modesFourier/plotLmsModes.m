directory = 'mur silvia\modesFourier\';

files{1} = 'LMS_ModalShapes_P0.xlsx';
files{2} = 'LMS_ModalShapes_P6.xlsx';
files{3} = 'LMS_ModalShapes_P7.xlsx';

P = [0 6 7];

for k = 1:3
    tab = xlsread([directory, files{k}]);
    
    p = P(k);
    
    for mode = 1:3
        L0 = 11*(mode-1);
        freq = tab(L0 + 1, 2);
        
        shape = tab((L0+2):(L0+10), 3) + 1i*tab((L0+2):(L0+10), 4);
        
        shape = transpose(shape);
        shape = shape / sqrt(shape * shape.');
        shape = sign(real(shape(1))) * shape;
        
        title = ['P', num2str(p), '_freq=', num2str(freq)];
        fig = plotModShape(real(shape), title);
        title2 = [title, '_complex'];
        fig2 = plotComplexModShape(shape, title2);
        
        % enregistrement
        savefig(fig, [directory, 'save\', title, '.fig']);
        saveas(fig, [directory, 'save\', title, '.png']);
        savefig(fig2, [directory, 'save\', title2, '.fig']);
        saveas(fig2, [directory, 'save\', title2, '.png']);
    end
end
