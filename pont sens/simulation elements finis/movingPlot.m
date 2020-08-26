function movingPlot(Xtot, t, L, essieux, c, pos_capteurs, timeCoeff)
% framerate interpolation
Tframes = (t(end) - t(1))/timeCoeff;
Nframes = length(t);

N = size(Xtot, 1);
dx = L / (N-1);
Xtot = Xtot / max(abs(Xtot), [], 'all');

fig = figure;
ax = axes(fig);
hold(ax, 'on');

xlim(ax, [-L/3, 4*L/3]);
ylim(ax, [-5, 5]);
xticks([0, L]);
xticklabels({'0', 'L'});
yticks([]);

    function [plt1, plt2, plt3, txt] = init_animation()
        cla(ax);
        
        plt1 = plot(ax, linspace(0, L, N), Xtot(:, 1), '-+');
        essieux_1 = essieux + c*t(1);
        plt2 = scatter(ax, essieux_1,...
            Xtot( max(min( floor(essieux_1/dx) + 1, N), 1), 1)' .* (1 - essieux_1/dx + floor(essieux_1/dx))...
            + Xtot( max(min( floor(essieux_1/dx) + 2, N), 1), 1)' .* (essieux_1/dx - floor(essieux_1/dx)), 'r');
        plt3 = scatter(ax, pos_capteurs,...
            Xtot( max(min( floor(pos_capteurs/dx) + 1, N), 1), 1)' .* (1 - pos_capteurs/dx + floor(pos_capteurs/dx))...
            + Xtot( max(min( floor(pos_capteurs/dx) + 2, N), 1), 1)' .* (pos_capteurs/dx - floor(pos_capteurs/dx)), 'dm');
        txt = text(ax, 0, 4, sprintf('t = %.2f', t(1)));
    end

    function animation()
        try
            [plt1, plt2, plt3, txt] = init_animation();
            %         pause(1);
            
            time0 = 24*3600*now;
            it = 1;
            while true
                set(plt1, 'YData', Xtot(:, it));
                essieux_it = essieux + c*t(it);
                set(plt2, 'XData', essieux_it, 'YData',...
                    Xtot( max(min( floor(essieux_it/dx) + 1, N), 1), it)' .* (1 - essieux_it/dx + floor(essieux_it/dx))...
                    + Xtot( max(min( floor(essieux_it/dx) + 2, N), 1), it)' .* (essieux_it/dx - floor(essieux_it/dx)));
                set(plt3, 'YData',...
                    Xtot( max(min( floor(pos_capteurs/dx) + 1, N), 1), it)' .* (1 - pos_capteurs/dx + floor(pos_capteurs/dx))...
                    + Xtot( max(min( floor(pos_capteurs/dx) + 2, N), 1), it)' .* (pos_capteurs/dx - floor(pos_capteurs/dx)));
                set(txt, 'String', sprintf('t = %.2f', t(it)));
                drawnow;
%                 pause(0.001);
                
                if it == Nframes
                    break
                end
                it = round((24*3600*now-time0)/Tframes * Nframes) + 1;
                if it > Nframes
                    it = Nframes;
                end
            end
        catch
        end
    end

init_animation();
animation

set(ax, 'ButtonDownFcn', @(~, ~) animation);

end

