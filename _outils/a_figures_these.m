%% moyen format

fig = gcf;
fig.Position(3:4) = [460 320];
set(fig, 'Renderer', 'painters');

%% grand format

fig = gcf;
fig.Position(3:4) = [560 420];
set(fig, 'Renderer', 'painters');

%% moyen petit (2 cote à cote)

fig = gcf;
fig.Position(3:4) = [350 250];
% set(fig, 'Renderer', 'painters');

%% petit (3 cote à cote)

fig = gcf;
fig.Position(3:4) = [260 210];
set(fig, 'Renderer', 'painters');


%%
%% xylabel normal

ax = gca;
ax.XLabel.FontSize = 11;
ax.XLabel.Interpreter = 'none';
ax.YLabel.FontSize = 11;
ax.YLabel.Interpreter = 'none';


%% xylabel latex

xlabel('$\xi$');
ylabel('$\psi(\xi)$');

ax = gca;
ax.XLabel.FontSize = 12;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.FontSize = 12;
ax.YLabel.Interpreter = 'latex';


%% mise en forme figures CWT

fig = gcf;
ax = gca;

set(fig, 'Renderer', 'opengl');

ax.XLabel.String = 'Temps [s]';
ax.YLabel.String = 'Fréquence [Hz]';

ax.Layer = 'top';
ax.TickDir = 'in';
ax.Box = 'on';

O = findobj(ax, '-not', 'Type', 'axes');
for k = 1:length(O)
    switch O(k).Type
        case 'surface'
            O(k).CData = O(k).CData;
            O(k).FaceColor = 'interp';
        case 'line'
            O(k).LineWidth = 1.5;
    end
    O(k).ZData = 0 * O(k).ZData;
end

ax.ZScale = 'linear';
ax.ZLim = [0 1];

%% enregistrement figures CWT

fig = gcf;
ax = gca;

% set(fig, 'Renderer', 'painters');

waitfor(msgbox('600 dpi => .png'));

L = findobj(ax, 'Type', 'line');
S = findobj(ax, '-not', 'Type', 'axes', '-not', 'Type', 'line');
Leg = findobj(fig, 'Type', 'legend');
set(Leg, 'AutoUpdate', 'off');
LongL = gobjects(0);
LongLXData = {};

% lignes avec bcp de donnees
k = 1;
while k <= length(L)
    if length(L(k).XData) > 1e5
        LongL(end+1) = L(k);
        L = [L(1:k-1), L(k+1:end)];
    else
        k = k+1;
    end
end

% enregistrement axes
for k = 1:length(S)
    S(k).Visible = 'off';
end
for k = 1:length(LongL)
    LongLXData{end+1} = LongL(k).XData;
    LongL(k).XData(:) = nan;
end
fig.Color = 'none';
ax.Color = 'none';
drawnow
try
    export_fig(gcf, '.eps', '-nocrop') % bug mais indispensable !
catch
end
waitfor(msgbox('.eps'));

% enregistrement axes 2
ax.XRuler.Axle.Visible = 'off';
ax.YRuler.Axle.Visible = 'off';
TickLength0 = ax.TickLength;
ax.TickLength = [0 0];
waitfor(msgbox('2.eps'));

% enregistrement image
ax.XRuler.Axle.Visible = 'on';
ax.YRuler.Axle.Visible = 'on';
ax.TickLength = TickLength0;
for k = 1:length(S)
    S(k).Visible = 'on';
end
for k = 1:length(LongL)
    LongL(k).XData = LongLXData{k};
end
for k = 1:length(L)
    L(k).Visible = 'off';
end
try
    Leg.Visible = 'off';
end
set(ax, 'Position', ax.Position); % pour fixer les axes
axXColor0 = ax.XColor;
axYColor0 = ax.YColor;
ax.XColor = 'none';
ax.YColor = 'none';
fig.Color = 'w';
waitfor(msgbox('2.png'));

% enregistrement image avec bords et ticks
ax.XColor = axXColor0;
ax.YColor = axYColor0;
ax.XTickLabel = {};
ax.YTickLabel = {};
xlabel('');
ylabel('');
waitfor(msgbox('4.png'));

% enregistrement image avec bords
xticks([]);
yticks([]);
waitfor(msgbox('3.png'));


%% reechantillonnage /10 lignes

n = 3;

fig = gcf;
ax = gca;

set(fig, 'Renderer', 'painters');

L = findobj(ax, 'Type', 'line');

for k = 1:length(L)
    set(L(k), 'XData', L(k).XData(1:n:end), 'YData', L(k).YData(1:n:end));
end














