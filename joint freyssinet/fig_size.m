
L = 1100;
H = 340;
% L = 1000;
% H = 500;
% L = 420;
% H = 380;
L = 600;
H = 340;
L = 700;
H = 550;

%%
fig = gcf;
ax = gca;

set(fig, 'Units', 'characters');
set(ax, 'Units', 'characters');
margins = [ax.Position([1 2]), fig.InnerPosition([3 4]) - ax.Position([1 2]) - ax.Position([3 4])];
set(fig, 'Units', 'pixels');
set(ax, 'Units', 'normalized');

%%
fig.Position(3) = L;
fig.Position(4) = H;

% set(fig, 'Units', 'characters');
% set(ax, 'Units', 'characters');
% ax.Position = [margins([1 2]), fig.InnerPosition([3 4]) - margins([1 2]) - margins([3 4])];
% set(fig, 'Units', 'pixels');
% set(ax, 'Units', 'normalized');