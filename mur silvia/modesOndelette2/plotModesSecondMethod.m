directory = 'mur silvia\modesOndelette2\';

load([directory, 'ModesWvlt2.mat']);

P = [0 6 7];

delete([directory, 'save\*']);

stringInfo = '';

for k = [1, 3]
    p = P(k);
    
    for mode = 1:size(ModesWvlt2, 2)
        if isempty(ModesWvlt2(k, mode).freq)
            continue
        end
        
        shape = ModesWvlt2(k, mode).shape;
        
        title = ['P', num2str(p), 'M', num2str(mode), 'CWT2'];
        fig = plotModShape(real(shape), title);
%         fig0 = plotModShape(imag(shape), title);
        title2 = [title, '_complex'];
        fig2 = plotComplexModShape(shape, title2);
        
        % enregistrement figures
        savefig(fig, [directory, 'save\', title, '.fig']);
        saveas(fig, [directory, 'save\', title, '.png']);
        saveas(fig, [directory, 'save\', title, '.eps'], 'epsc');
        savefig(fig2, [directory, 'save\', title2, '.fig']);
        saveas(fig2, [directory, 'save\', title2, '.png']);
        saveas(fig2, [directory, 'save\', title2, '.eps'], 'epsc');
    end
    
    stringInfo = [stringInfo, 'P', num2str(p), newline];
    stringInfo = [stringInfo, ModesWvlt2(k, 2).choixQ, newline];
    stringInfo = [stringInfo, ModesWvlt2(k, 2).choixNoise, newline, newline];
end

infoFile = fopen([directory, 'save\info.txt'], 'w');
fwrite(infoFile, stringInfo);
fclose(infoFile);




