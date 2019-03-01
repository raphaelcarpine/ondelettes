function animate(objet, t, x)
%animate Summary of this function goes here
%   Detailed explanation goes here
t = t-t(1);
N = length(t);

deltaT = 15;
dT = 0.01;
gamma = t(end)/deltaT;

x = x(:, 1:int16(size(x,2)/2));

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
for k=1:n
    if isa(objet, 'Structure') && k<=objet.ddl
        alpha = 1/max(abs(X(:, k)));
        draw = plot(ax, k-1, 0, 'o');
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
            draw2 = plot(ax, ddl-1, -L, 'o');
            updatek = @(Xt) set(draw1, 'XData', [ddl-1+0.5*struc*alpha*Xt(ddl) ddl-1+0.5*struc*alpha*Xt(ddl)+L*sin(Xt(k))], 'YData', [0 -L*cos(Xt(k))]);
            updatek1 = @(Xt) set(draw2, 'XData', ddl-1+0.5*struc*alpha*Xt(ddl)+L*sin(Xt(k)), 'YData', -L*cos(Xt(k)));
            updates{end+1} = updatek;
            updates{end+1} = updatek1;
        elseif isa(tmd, 'TMDmasseressort')
            alpha = 1/max(abs(X(:, ddl)));
            draw1 = plot(ax, [ddl-1 ddl-1 ddl-1 ddl-1], [0 -0.5 -0.5 -1]);
            ax.ColorOrderIndex = ax.ColorOrderIndex-1;
            draw2 = plot(ax, ddl-1, -1, 'o');
            updatek = @(Xt) set(draw1, 'XData', [ddl-1+0.5*struc*alpha*Xt(ddl), ddl-1+0.5*struc*alpha*Xt(ddl), ddl-1+0.5*struc*alpha*Xt(ddl)+0.5*alpha*Xt(k), ddl-1+0.5*struc*alpha*Xt(ddl)+0.5*alpha*Xt(k)]);
            updatek1 = @(Xt) set(draw2, 'XData', ddl-1+0.5*struc*alpha*Xt(ddl)+0.5*alpha*Xt(k));
            updates{end+1} = updatek;
            updates{end+1} = updatek1;
        end
    end
end

    function show()
        try
            for Xt=X.'
                for kupdate = 1:length(updates)
                    update = updates{kupdate};
                    update(Xt);
                end
                drawnow limitrate;
                pause(dT);
            end
        catch
        end
    end

drawnow;
% show();
end

