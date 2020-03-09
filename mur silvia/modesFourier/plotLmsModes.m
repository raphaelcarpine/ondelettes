directory = 'mur silvia\modesFourier\';

load([directory, 'ModesLMS.mat']);

P = [0 6 7];

delete([directory, 'save\*']);

for k = 1:3
    tab = xlsread([directory, files{k}]);
    
    p = P(k);
    
    mode = 1;
    while mode <= size(ModesLMS, 2) && ~isempty(ModesLMS(k, mode).freq)
        freq = ModesLMS(k, mode).freq;
        shape = ModesLMS(k, mode).shape;
        damping = ModesLMS(k, mode).damping;
        
        title = ['P', num2str(p), 'M', num2str(mode), '_freq=', num2str(freq), '_damping=', num2str(100*damping)];
        fig = plotModShape(real(shape), title);
%         fig0 = plotModShape(imag(shape), title);
        title2 = [title, '_complex'];
        fig2 = plotComplexModShape(shape, title2);
        
        % enregistrement figures
        savefig(fig, [directory, 'save\', title, '.fig']);
        saveas(fig, [directory, 'save\', title, '.png']);
        saveas(fig, [directory, 'save\', title, '.eps']);
        savefig(fig2, [directory, 'save\', title2, '.fig']);
        saveas(fig2, [directory, 'save\', title2, '.png']);
        saveas(fig2, [directory, 'save\', title2, '.eps']);
        
        mode = mode + 1;
    end
end