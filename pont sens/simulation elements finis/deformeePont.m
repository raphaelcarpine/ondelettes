function fctDeformee = deformeePont(L, pos_capteurs)
%DEFORMEEPONT Summary of this function goes here
%   Detailed explanation goes here

n = length(pos_capteurs);

pos_capteurs = [0, pos_capteurs, L];

% préparation points abcsisse
frac_pos_capteurs = pos_capteurs;
frac_pos_capteurs = frac_pos_capteurs/L;
frac_pos_capteurs = rats(frac_pos_capteurs);
frac_pos_capteurs = strsplit(frac_pos_capteurs, ' ');
ifrac = 1;
while ifrac <= length(frac_pos_capteurs)
    if isempty(frac_pos_capteurs{ifrac})
        frac_pos_capteurs = {frac_pos_capteurs{1:ifrac-1}, frac_pos_capteurs{ifrac+1:end}};
    else
        ifrac = ifrac + 1;
    end
end
for ifrac = 1:length(frac_pos_capteurs)
    frac_pos_capteurs{ifrac} = strsplit(frac_pos_capteurs{ifrac}, '/');
    if strcmp(frac_pos_capteurs{ifrac}{1}, '1')
        frac_pos_capteurs{ifrac}{1} = 'L';
    elseif ~strcmp(frac_pos_capteurs{ifrac}{1}, '0')
        frac_pos_capteurs{ifrac}{1} = [frac_pos_capteurs{ifrac}{1}, 'L'];
    end
    if length(frac_pos_capteurs{ifrac}) > 1
        frac_pos_capteurs{ifrac} = {frac_pos_capteurs{ifrac}{1}, '/', frac_pos_capteurs{ifrac}{2}};
    end
    frac_pos_capteurs{ifrac} = [frac_pos_capteurs{ifrac}{:}];
end



    function fctDeformeePont(shape, figTitle)
        if nargin < 2
            figTitle = '';
        end
        
        if length(shape) ~= n
            error('wrong shape array size');
        end
        
        if iscolumn(shape)
            shape = shape.';
        end
        shape = [0, shape, 0];
        
        fig = figure('Name', figTitle);
        ax = axes(fig);
        hold(ax, 'on');
        plot(ax, pos_capteurs, zeros(size(pos_capteurs)), '--o', 'Color', 0.5*[1 1 1]);
        plot(ax, pos_capteurs, shape, '-o');
        xlim(ax, L*[-0.1, 1.1]);
        ylim(ax, [min(shape), max(shape)] + 0.1*(max(shape)-min(shape))*[-1 1]);
        xlabel(ax, 'x');
        ylabel(ax, '\phi');
        xticks(ax, pos_capteurs);
        xticklabels(ax, frac_pos_capteurs);
    end

fctDeformee = @fctDeformeePont;

end

