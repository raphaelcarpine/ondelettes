function [fctDefModale, fctDefModaleAnimation] = defModalePontLin(L, pos_capteurs, nameFig)
%DEFMODALES Summary of this function goes here
%   pile : 1-10
%   chNames = {'acc1x', 'acc1y', 'acc1z', 'acc2x', 'acc2y', 'acc3x',... }

if nargin < 3
    nameFig = '';
end

%% fonction deformee

    function fig = fctDefModale0(animated, deformee, nameFig2)
        if nargin < 3
            nameFig2 = nameFig;
        end
        
        if length(deformee) ~= length(pos_capteurs)
            error('length(deformee) ~= length(pos_capteurs)');
        end
        if ~animated && any(imag(deformee) ~= 0)
            error('complex mode shape');
        end
        
        deformee = deformee * sign(max(real(deformee)));
        deformee = deformee / max(abs(deformee));
        
        %%
        
        fig = figure('Name', nameFig2);
        fig.Position(3:4) = [560 200];
        fig.UserData = {L, pos_capteurs, deformee}; % save
        fig.Position(3:4) = [250 350];
        ax = axes(fig);
        ax.Position = [0 0 1 1];
        hold(ax, 'on');
        
        % pont
        col0 = 0.5*[1 1 1];
        plot(ax, [0 L], [0 0], 'Color', col0);
        
        % deformée
        def0 = real(deformee);
        col1 = 0.*[1 1 1];
        defPlot = plot(ax, pos_capteurs, def0, '-o', 'Color', col1);
        
        % display
        set(ax, 'XLim', L*([0 1] + 0.1*[-1 1]));
        set(ax, 'YLim', 2*[-1 1]);
        set(ax,'visible','off');
        
        %% animation
        if ~animated
            return
        end
        
        freq = 0.6; % display frequency
        fps = 25; % nb of frames
        
        function animation
            Nframes = ceil(fps/freq);
            tau = 1/fps;
            t = (0:Nframes-1)/Nframes;
            tic;
            while fig == get(groot, 'CurrentFigure')
                try
                    k = floor(toc/tau) +1;
                    realShape = real(deformee * exp(2i*pi*t(mod(k-1, Nframes) +1)));
                    set(defPlot, 'YData', realShape);
                    drawnow
                    pause((k+1)*tau - toc);
                catch
                    break
                end
            end
            
            % end
            try
                realShape = real(deformee);
                set(defPlot, 'YData', realShape);
                drawnow
            catch
            end
        end
        
        clickTxt = annotation(fig, 'textbox', [0 0 1 1], 'String', 'Click to start animation', 'FitBoxToText', 'off',...
            'FontSize', 16, 'BackgroundColor', 'white', 'FaceAlpha', 0.5, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        set(get(ax,'Children'),'HitTest','off')
        
        
        function clickFcn(~,~)
            delete(clickTxt);
            animation;
        end
        
        fig.WindowButtonDownFcn = @clickFcn;
        fig.Interruptible = 'off';
        
    end


fctDefModale = @(varargin) fctDefModale0(false, varargin{:});
fctDefModaleAnimation = @(varargin) fctDefModale0(true, varargin{:});


end

