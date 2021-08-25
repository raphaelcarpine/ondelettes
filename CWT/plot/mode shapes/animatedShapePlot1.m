function animatedShapePlot1(shape, figTitle)

freq = 0.7; % display frequency
fps = 25; % nb of frames

if nargin < 1
    warning('shape mising');
    shape = randn(1, 5) + 0.3i*rand(1, 5);
end

if nargin < 2
    figTitle = '';
end

%% plot

n = length(shape);

fig = figure('Name', figTitle);
ax = axes(fig);
plt = stem(ax, real(shape));
xlim(ax, [0, n+1]);
ylim(ax, max(abs(shape)) * [-1.2, 1.2]);
xlabel(ax, 'dof');
ylabel(ax, '\phi');
xticks(ax, 1:n);

%% animation

    function animation
        Nframes = ceil(fps/freq);
        tau = 1/fps;
        t = (0:Nframes-1)/Nframes;
        k = 1;
        tic;
        while fig == get(groot,'CurrentFigure')
            try
                k = floor(toc/tau) +1;
                realShape = real(shape * exp(2i*pi*t(mod(k-1, Nframes) +1)));
                set(plt, 'YData', realShape);
                drawnow
                pause((k+1)*tau - toc);
            catch
                break
            end
        end
        
        % end
        try
            set(plt, 'YData', real(shape));
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

