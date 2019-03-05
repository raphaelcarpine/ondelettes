function animate(objet, t, x)
%animate Summary of this function goes here
%   Detailed explanation goes here
t = t-t(1);
N = length(t);

deltaT = 15;
dT = 0.01;
gamma = t(end)/deltaT;


X = zeros(0);

n = 1;
for T=0:dT:deltaT
    while n<N && t(n+1)<=gamma*T
        n = n+1;
    end
    if n>=N
        break
    end
    c = (gamma*T-t(n))/(t(n+1)-t(n));
    X(end+1,:) = x(n,:)*c+x(n+1,:)*(1-c);
end


if isa(objet, 'Structure') && size(x, 2)==objet.ddl
    lignes = 1;
else
    lignes = 2;
end
if isa(objet, 'Structure')
    colonnes = objet.ddl;
else
    colonnes = 1;
end


fig = figure;
axis off;
% axis equal;
ax = gca;
set(ax,'DataAspectRatio',[1 1 1])
ax.Position = ax.OuterPosition;
hold(ax, 'on');
xlim([-1.5, colonnes+0.5]);
ylim([0.5-lignes, 0.5]);
set(fig, 'WindowButtonDownFcn', @(~, ~) show());

updates = {};
n = size(X, 2);
for k = 1:n
    if isa(objet, 'Structure') && k<=objet.ddl
        alpha = 1/max(abs(X(:, k)));
        draw = plot(ax, k-1, 0, 'o', 'MarkerSize',20);
        set(draw, 'MarkerFaceColor', get(draw, 'Color'));
        updatek = @(Xt) set(draw, 'XData', k-1+0.5*alpha*Xt(k));
        updates{end+1} = updatek;
    else
        if isa(objet, 'Structure')
            tmd = objet.TMDs{k-objet.ddl}{1};
            ddl = objet.TMDs{k-objet.ddl}{2};
            struc = 1;
        else
            tmd = objet;
            ddl = 1;
            struc = 0;
        end
        
        if isa(tmd, 'TMDpendule')
            L = 0.5*tmd.l/max(abs(X(:, ddl)));
            alpha = 1/max(abs(X(:, ddl)));
            draw1 = plot(ax, [ddl-1 ddl-1], [0 -L]);
            ax.ColorOrderIndex = ax.ColorOrderIndex-1;
            draw2 = plot(ax, ddl-1, -L, 'o', 'MarkerSize', 10);
            updatek = @(Xt) set(draw1, 'XData', [ddl-1+0.5*struc*alpha*Xt(ddl) ddl-1+0.5*struc*alpha*Xt(ddl)+L*sin(Xt(k))], 'YData', [0 -L*cos(Xt(k))]);
            updatek1 = @(Xt) set(draw2, 'XData', ddl-1+0.5*struc*alpha*Xt(ddl)+L*sin(Xt(k)), 'YData', -L*cos(Xt(k)));
            updates{end+1} = updatek;
            updates{end+1} = updatek1;
        elseif isa(tmd, 'TMDmasseressort')
            alpha = 1/max(abs(X(:, ddl)));
            draw1 = plot(ax, ressort(ddl-1, ddl-1, 0, -1, 1), ressort(ddl-1, ddl-1, 0, -1, 2));
            ax.ColorOrderIndex = ax.ColorOrderIndex-1;
            draw2 = plot(ax, ddl-1, -1, 'o', 'MarkerSize', 10);
            updatek = @(Xt) set(draw1, 'XData', ressort(ddl-1+0.5*struc*alpha*Xt(ddl), ddl-1+0.5*struc*alpha*Xt(ddl)+0.5*alpha*Xt(k), 0, -1, 1));
            updatek1 = @(Xt) set(draw2, 'XData', ddl-1+0.5*struc*alpha*Xt(ddl)+0.5*alpha*Xt(k));
            updates{end+1} = updatek;
            updates{end+1} = updatek1;
        end
    end
end

    function Update(updates, Xt)
        for kupdate = 1:length(updates)
            update = updates{kupdate};
            update(Xt);
        end
        drawnow limitrate;
    end

    function show()
        try
            for Xt=X.'
                Update(updates, Xt);
                pause(dT);
            end
        catch
        end
    end

% uistack(updates);
Update(updates, X(1, :));
% show();
end



function XY = ressort(x1, x2, y1, y2, dim)
alpha = 0.1;
n = 6;

L = x2-x1;
H = y2-y1;
XY = zeros(2, 4+2*n);

XY(:, 1) = [x1 y1];
XY(:, 2) = [x1 (y1+y2)/2];
XY(:, end-1) = [x2 (y1+y2)/2];
XY(:, end) = [x2 y2];

for k=1:n
    XY(1, 2*k+1) = x1 + L*(4*k-3)/(4*n);
    XY(1, 2*k+2) = x1 + L*(4*k-1)/(4*n);
    XY(2, 2*k+1) = (y1+y2)/2 + alpha*H;
    XY(2, 2*k+2) = (y1+y2)/2 - alpha*H;
end

XY = XY(dim, :);

end

