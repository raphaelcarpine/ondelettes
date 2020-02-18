h = 0.5;
M = -[0, 0, 1];

mu0 = 1;

l = 1;
L = 1;
N = 1000;
x = linspace(-L, L, N);
y = linspace(-l, l, N);


B0 = @(r) mu0/(4*pi) * ( 3*(M*r.')*r - M ) / (r*r.')^(3/2);

B0z = nan(N, N);
for ix = 1:N
    for iy = 1:N
        B0z(iy, ix) = B0([x(ix), y(iy), -h]) * [0; 0; 1];
    end
end

fig = figure;
ax = axes(fig);
plt = pcolor(ax, x, y, B0z);
axis equal
axis tight
set(ax,'xcolor','none')
ax.XAxis.Label.Color=[0 0 0];
ax.XAxis.Label.Visible='on';
set(ax,'ycolor','none')
ax.YAxis.Label.Color=[0 0 0];
ax.YAxis.Label.Visible='on';
xlabel(ax, 'x')
ylabel(ax, 'y')
shading flat
colormap(ax, jet);
colorbar(ax);

Xlim = 3;
deltaX = 0.05;
for x0 = -Xlim:deltaX:Xlim
    B0z = nan(N, N);
    for ix = 1:N
        for iy = 1:N
            B0z(iy, ix) = B0([x(ix) - x0, y(iy), -h]) * [0; 0; 1];
        end
    end
%     set(plt, 'ZData', B0z);
%     drawnow
    plt = pcolor(ax, x, y, B0z);
    axis equal
    axis tight
    set(ax,'xcolor','none')
    ax.XAxis.Label.Color=[0 0 0];
    ax.XAxis.Label.Visible='on';
    set(ax,'ycolor','none')
    ax.YAxis.Label.Color=[0 0 0];
    ax.YAxis.Label.Visible='on';
    xlabel(ax, 'x')
    ylabel(ax, 'y')
    shading flat
    colormap(ax, jet);
    colorbar(ax);
    caxis(ax, [0, 0.16]);
    drawnow
    if x0 == -Xlim
        gif('maquetteTMD/courbes cours/champB.gif');
    else
        gif;
    end
end



%%

fig = figure;
ax = axes(fig);

for x0 = -Xlim:deltaX:Xlim
    B0z = nan(N, N);
    for ix = 1:N
        for iy = 1:N
            B0z(iy, ix) = B0([x(ix) - x0, y(iy), -h]) * [0; 0; 1];
        end
    end
    dB0z = diff(B0z')';
    dB0z = [dB0z(:,1), dB0z] + [dB0z, dB0z(:,end)];
    
    plt = pcolor(ax, x, y, - dB0z);
    axis equal
    axis tight
    set(ax,'xcolor','none')
    ax.XAxis.Label.Color=[0 0 0];
    ax.XAxis.Label.Visible='on';
    set(ax,'ycolor','none')
    ax.YAxis.Label.Color=[0 0 0];
    ax.YAxis.Label.Visible='on';
    xlabel(ax, 'x')
    ylabel(ax, 'y')
    shading flat
    colormap(ax, jet);
    colorbar(ax);
    caxis(ax, [-1.1, 1.1]*1e-3);
    drawnow
    
    if x0 == -Xlim
        gif('maquetteTMD/courbes cours/rotE.gif');
    else
        gif;
    end
end
