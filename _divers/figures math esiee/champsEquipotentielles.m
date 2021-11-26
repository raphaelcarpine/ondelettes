x = linspace(-3, 3, 100);
y = linspace(-3, 3, 100);

[X, Y] = meshgrid(x, y);

F = X.^2 + Y.^2; Ftitle = 'F(x,y) = x² + y²';
% F = X.^2 + (Y-1).^2; Ftitle = 'F(x,y) = x² + (y-1)²';
% F = X; Ftitle = 'F(x,y) = x';
% F = Y.^2 + 2*Y; Ftitle = 'F(x,y) = y² + 2y';
% F = Y./(1+X.^2); Ftitle = 'F(x,y) = y/(1+x^2)';

% champ

figure;
surf(X, Y, F, 'EdgeColor', 'none');
hold on
title(Ftitle);
xlabel('x');
ylabel('y');
zlabel('F(x,y)');
set(gca, 'DataAspectRatio', [diff(get(gca, 'XLim')) diff(get(gca, 'XLim')) diff(get(gca, 'ZLim'))]);
view(0, 90);

colormap jet
hcb = colorbar;
hcb.Title.String = 'F(x,y)';

equi = [];
equiSurf = [];

%% equipotentielle
delete(equi);
delete(equiSurf);

lambda = 5;

if isnan(lambda)
    return
end


% equiSurf = surf(X, Y, lambda*ones(size(F)), 'EdgeColor', 'none', 'FaceColor', [0 0 0], 'FaceAlpha', 0.3);


% lignes intersection
z1 = F;
z2 = lambda*ones(size(F));
zdiff = z1 - z2;
C = contours(x, y, zdiff, [0 0]);
% Extract the x- and y-locations from the contour matrix C.
xL = C(1, 2:end);
yL = C(2, 2:end);
% Interpolate on the first surface to find z-locations for the intersection
% line.
zL = interp2(x, y, z1, xL, yL);
% Visualize the line.
equi = line(xL, yL, zL, 'Color', 'k', 'LineWidth', 3, 'LineStyle', '--');

% lambda sur la colorbar
set(hcb, 'XTickMode', 'auto', 'XTickLabelMode', 'auto');
XTickCB = get(hcb, 'XTick');
XTicklabelCB = get(hcb, 'XTickLabel');
XTicklabelCB = XTicklabelCB(XTickCB ~= lambda);
XTickCB = XTickCB(XTickCB ~= lambda);
XTickCB = [XTickCB, lambda];
XTicklabelCB{end+1} = '\lambda';
[XTickCB, Ixtick] = sort(XTickCB);
XTicklabelCB = XTicklabelCB(Ixtick);
set(hcb, 'XTick', XTickCB);
set(hcb, 'XTickLabel', XTicklabelCB);