function movingPlotVehiclesGif(Xtot, T, L, t_vehicles_left, t_vehicles_right, m_vehicles_left, m_vehicles_right,...
    c_vehicles_left, c_vehicles_right, pos_capteurs, pos_nonlin, nonlin_reached, timeCoeff)

fileName = 'NL_GLOB.gif';
frameInterp = 1;
ti = 10;
tf = 100;
plotTime = 0;

plotCapt = false;

if isempty(nonlin_reached) % lin
    nonlin_reached = nan(size(T));
elseif size(nonlin_reached, 1) > 1 % nonlin global
    pos_nonlin = linspace(0, L, size(Xtot, 1));
    nonlin_reached = [false(1, length(T)); nonlin_reached; false(1, length(T))]; % ddl manquants
end
nonlin_reached_nans = nonlin_reached ./ nonlin_reached - 1; % 0 si true, nan sinon

% framerate interpolation
Tframes = (T(end) - T(1))/timeCoeff;
Nframes = length(T);

N = size(Xtot, 1);
dx = L / (N-1);
Xtot = 1e3 * Xtot; % [mm]

fig = figure;
fig.Renderer = 'painters';
fig.Color = 'white';
fig.Position(3:4) = [560 200];
ax = axes(fig);
ax.Units = 'pixels';
ax.Position = [1 1 560 200] + 40*[1 0 -1 0] + 40*[0 1 0 -1] + 20*[0 0 -1 -1];
hold(ax, 'on');

xlim(ax, L*[0 1] + 5*[-1 1]);
ylimmax = 1/2*ceil(2*1.2 * max(abs(Xtot), [], 'all'));
ylim(ax, ylimmax*[-1 1]);
% xticks([0, L]);
% xticklabels({'0', 'L'});
% yticks([]);
% daspect(ax, [1 1 1]);
xlabel('x [m]');
ylabel('w [mm]');
ax.Color = 'none';
% axis tight
ax.Clipping = 'off';
% axis off

ind_nonlin = round(pos_nonlin/dx) + 1;

    function [plt1, plt2, plt3l, plt3r, plt4, txt] = init_animation()
        cla(ax); % clear axes
        
        plt1 = plot(ax, linspace(0, L, N), Xtot(:, 1), 'k-+');
        vehicles_1l = c_vehicles_left .* (T(1) - t_vehicles_left);
        vehicles_1r = L - c_vehicles_right .* (T(1) - t_vehicles_right);
        if plotCapt
            plt2 = scatter(ax, pos_capteurs,...
                Xtot( max(min( floor(pos_capteurs/dx) + 1, N), 1), 1)' .* (1 - pos_capteurs/dx + floor(pos_capteurs/dx))...
                + Xtot( max(min( floor(pos_capteurs/dx) + 2, N), 1), 1)' .* (pos_capteurs/dx - floor(pos_capteurs/dx)), 'dm');
            
        else
            plt2 = [];
        end
        plt3l = scatter(ax, vehicles_1l,...
            Xtot( max(min( floor(vehicles_1l/dx) + 1, N), 1), 1)' .* (1 - vehicles_1l/dx + floor(vehicles_1l/dx))...
            + Xtot( max(min( floor(vehicles_1l/dx) + 2, N), 1), 1)' .* (vehicles_1l/dx - floor(vehicles_1l/dx)),...
            max(100 * (m_vehicles_left/mean(m_vehicles_left)), 1), 'b>', 'MarkerFaceColor', 'flat');
        plt3r = scatter(ax, vehicles_1r,...
            Xtot( max(min( floor(vehicles_1r/dx) + 1, N), 1), 1)' .* (1 - vehicles_1r/dx + floor(vehicles_1r/dx))...
            + Xtot( max(min( floor(vehicles_1r/dx) + 2, N), 1), 1)' .* (vehicles_1r/dx - floor(vehicles_1r/dx)),...
            max(100 * (m_vehicles_right/mean(m_vehicles_left)), 1), 'b<', 'MarkerFaceColor', 'flat');
        plt4 = scatter(ax, pos_nonlin, Xtot(ind_nonlin, 1) + nonlin_reached_nans(:, 1), 300, '.',...
            'MarkerEdgeColor', 'red', 'LineWidth', 0.5);
        txt = text(ax, 0, 0.6*ylimmax, sprintf('t = %.2f', T(1)));
        if ~plotTime
            txt.Visible = 'off';
        end
    end

    function update_fig(kt)
        set(plt1, 'YData', Xtot(:, kt));
        vehicles_itl = c_vehicles_left .* (T(kt) - t_vehicles_left);
        vehicles_itr = L - c_vehicles_right .* (T(kt) - t_vehicles_right);
        if plotCapt
            set(plt2, 'YData',...
                Xtot( max(min( floor(pos_capteurs/dx) + 1, N), 1), kt)' .* (1 - pos_capteurs/dx + floor(pos_capteurs/dx))...
                + Xtot( max(min( floor(pos_capteurs/dx) + 2, N), 1), kt)' .* (pos_capteurs/dx - floor(pos_capteurs/dx)));
        end
        set(plt3l, 'XData', vehicles_itl, 'YData',...
            Xtot( max(min( floor(vehicles_itl/dx) + 1, N), 1), kt)' .* (1 - vehicles_itl/dx + floor(vehicles_itl/dx))...
            + Xtot( max(min( floor(vehicles_itl/dx) + 2, N), 1), kt)' .* (vehicles_itl/dx - floor(vehicles_itl/dx)));
        set(plt3r, 'XData', vehicles_itr, 'YData',...
            Xtot( max(min( floor(vehicles_itr/dx) + 1, N), 1), kt)' .* (1 - vehicles_itr/dx + floor(vehicles_itr/dx))...
            + Xtot( max(min( floor(vehicles_itr/dx) + 2, N), 1), kt)' .* (vehicles_itr/dx - floor(vehicles_itr/dx)));
        set(plt4, 'YData', Xtot(ind_nonlin, kt) +  + nonlin_reached_nans(:, kt));
        set(txt, 'String', sprintf('t = %.2f', T(kt)));
        drawnow;
    end

    function animation(~, ~)
        try
            gif(fileName, 'frame', gcf, 'DelayTime', (T(end)-T(1))/(length(T)-1)*frameInterp);
            
            kti = 1;
            while kti < length(T) && T(kti) < ti
                kti = kti + 1;
            end
            ktf = 0;
            while ktf < length(T) && T(ktf+1) <= tf
                ktf = ktf + 1;
            end
            for kt = kti:frameInterp:ktf
                update_fig(kt);
                gif;
            end
        catch
        end
    end

[plt1, plt2, plt3l, plt3r, plt4, txt] = init_animation();

animation;

end

