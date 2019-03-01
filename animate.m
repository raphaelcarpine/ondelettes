function animate(objet, t, x)
%animate Summary of this function goes here
%   Detailed explanation goes here
t = t-t(1);
N = length(t);

deltaT = 15;
dT = 0.04;
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
    X(end+1,:) = 0.8*(x(n,:)*c+x(n+1,:)*(1-c));
end


if isa(objet, 'Structure') && size(x, 1)==objet.ddl
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
xlim([-1 colonnes]);
ylim([0.5-lignes 0.5]);
set(fig, 'WindowButtonDownFcn', @(~, ~) show());

updates = {};
n = size(X, 2);
X(1, end+1) = 0;
for k=1:n
    if isa(objet, 'Structure') && k<=objet.ddl
        X(:, k) = X(:, k)/max(abs(X(:, k)));
        draw = plot(ax, k-1, 0, 'o');
        updatek = @(Xt) set(draw, 'XData', k-1+0.4*Xt(k));
        updates{end+1} = updatek;
    else
        if isa(objet, 'Structure')
            tmd = objet.TMDs{k-objet.ddl}{1};
            ddl = objet.TMDs{k-objet.ddl}{2};
        else
            tmd = objet;
            ddl = n+1;
        end
        
        if isa(tmd, 'TMDpendule')
            draw1 = plot(ax, [ddl-1 ddl-1], [0 -1]);
%             draw2 = plot(ax, ddl-1, -1, 'o');
            updatek = @(Xt) set(draw1, 'XData', [ddl-1+0.4*Xt(ddl) ddl-1+0.4*Xt(ddl)+sin(Xt(k))], 'YData', [0 -cos(Xt(k))]);
%             set(draw2, 'XData', ddl-1+0.4*Xt(ddl)+0.4*sin(Xt(k)), 'YData', -0.4*cos(Xt(k)));
            updates{end+1} = updatek;
        end
    end
end

    function show()
        for Xt=X.'
            for kupdate = 1:length(updates)
                update = updates{kupdate};
                update(Xt);
            end
            drawnow limitrate;
            pause(dT);
        end
    end

show();
end

