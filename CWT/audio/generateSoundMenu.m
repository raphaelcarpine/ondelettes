function generateSoundMenu()
%GENERATESOUNDMENU Summary of this function goes here
%   Detailed explanation goes here
Fs = 48000;
T = 1;
T0 = 0.5;

f = 200;
Q = 1000;


%%
N = round(T*Fs);
N0 = round(T0*Fs);
T = N/Fs;

%% display

L = 100;
H = 20;

fig = figure('Name', 'Sound generator', 'numbertitle', 'off');
fig.Units = 'characters';
fig.Position(3) = L;
fig.Position(4) = H;
fig.MenuBar = 'none';
fig.ToolBar = 'none';
% fig.Resize = false;

ax = axes(fig);
ax.Units = 'characters';
ax.Position = [2, 2, L-30, H-3];
set(ax, 'YTick', []);

fftPlt = plot(ax, nan, nan);
updateFftPlt();

uicontrol(fig, 'Style', 'text', 'String', 'f = ', 'HorizontalAlignment', 'right',...
    'Units', 'characters', 'Position', [L-24, H - 4, 5, 1.5]);
fInput = uicontrol(fig, 'Style', 'edit', 'String', f,...
    'Units', 'characters', 'Position', [L-19, H - 4, 7, 1.5]);
uicontrol(fig, 'Style', 'text', 'String', 'Q = ', 'HorizontalAlignment', 'right',...
    'Units', 'characters', 'Position', [L-24, H - 6, 5, 1.5]);
QInput = uicontrol(fig, 'Style', 'edit', 'String', Q,...
    'Units', 'characters', 'Position', [L-19, H - 6, 7, 1.5]);

butt = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Generate sound',...
    'Units', 'characters', 'Position', [L-25, H - 12, 20, 2]);


fInput.Callback = @updateInputs;
QInput.Callback = @updateInputs;
butt.Callback = @buttonCallback;

%% callbacks

    function updateInputs(~,~)
        f2 = str2double(get(fInput, 'String'));
        Q2 = str2double(get(QInput, 'String'));
        
        if isnan(f2)
            warning('incorrect input');
        else
            f = f2;
        end
        if isnan(Q2)
            warning('incorrect input');
        else
            Q = Q2;
        end
        
        updateFftPlt()
    end

    function updateFftPlt()
        zeta = 1/(2*Q);
        fn = f / sqrt(1-2*zeta^2);
        
        freqs = linspace(0, 2*f, 1000);
        Fx = 1./(1 - freqs.^2/fn^2 + 2i*zeta*freqs/fn);
        set(fftPlt, 'XData', freqs, 'YData', abs(Fx));
        set(ax, 'YTick', []);
        drawnow
    end





    function buttonCallback(~,~)
        updateInputs;
        x = generateSound();
        sound(x, Fs);
    end



    function x = generateSound()
        zeta = 1/(2*Q);
        fn = f / sqrt(1-2*zeta^2);
        
        freqs = (1/T) * [0:floor(N/2)-1, -N+floor(N/2):-1];
        Fx = 1./(1 - freqs.^2/fn^2 + 2i*zeta*freqs/fn);
        Fx = abs(Fx) .* exp(2i*pi*rand(size(Fx)));
        x = ifft(Fx);
        x = [(0:N0-1)/N0 * x(1), x, (N0-1:-1:0)/N0 * x(end)];
        x = real(x);
        x = x / max(abs(x));
    end

end

