directory = 'mur silvia\modesFourier\';

load([directory, 'ModesLMS.mat']);

P = [0 6 7];


for k = 1:3
    tab = xlsread([directory, files{k}]);
    
    p = P(k);
    
    for mode = 1:3
        freq = ModesLMS(k, mode).freq;
        shape = ModesLMS(k, mode).shape;
        
        title = ['P', num2str(p), '_freq=', num2str(freq)];
        fig = plotModShape(real(shape), title);
%         fig0 = plotModShape(imag(shape), title);
        title2 = [title, '_complex'];
        fig2 = plotComplexModShape(shape, title2);
        
        % enregistrement figures
        savefig(fig, [directory, 'save\', title, '.fig']);
        saveas(fig, [directory, 'save\', title, '.png']);
        savefig(fig2, [directory, 'save\', title2, '.fig']);
        saveas(fig2, [directory, 'save\', title2, '.png']);
    end
end