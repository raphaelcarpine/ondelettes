%% complexes Xn Xn+k

z0 = 0.6*exp(0.9i);
zeta = 0.01;
r0 = exp(-2*pi*zeta/sqrt(1-zeta^2));
z = r0*z0;
sigma = abs(z0) * 2*sqrt(pi*zeta/sqrt(1-zeta^2));

x = linspace(-0.2, 0.8, 300);
y = linspace(-0.1, 0.8, 300);

[X, Y] = meshgrid(x, y);


P = 1/(2*pi*sigma^2) * exp(-((X-real(z)).^2 + (Y-imag(z)).^2)/(2*sigma^2));

fig = figure;
% p = pcolor(X, Y, -P);
p = contourf(X, Y, P, 6);
% p = surf(X, Y, P);
% p.FaceColor = 'interp';
% p.EdgeColor = 'none';
colormap(flipud(hot)) %bone gray hot autumn
colorbar
% colorbar('XTick', 0)
caxis([0 inf])

yline(0, '-', 'Re');
xline(0, '-', 'Im', 'LabelHorizontalAlignment', 'left', 'LabelOrientation', 'horizontal');

color1 = [0.3 0.6 1];
color2 = [0.2 0.9 0];

hold on
plot(real(z0), imag(z0), '+', 'Color', color1, 'LineWidth', 1.5, 'MarkerSize', 6);
plot(real(z), imag(z), '+', 'Color', color2, 'LineWidth', 1.5, 'MarkerSize', 6);

% text(real(z0), imag(z0), {'  $\mathbf{X[n]}$', ''}, 'Color', color1, 'Interpreter', 'latex', 'FontSize', 12);
% text(real(z), imag(z), {'', '', '$\mathbf{r_0 X[n]}$ '}, 'Color', color2,...
%     'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'FontSize', 12, 'FontName', 'arial');
text(real(z0), imag(z0), {'  X[n]', ''}, 'Color', color1,...
    'Interpreter', 'tex', 'FontSize', 11, 'FontWeight', 'bold');
text(real(z), imag(z), {'', 'r_0X[n] '}, 'Color', color2, 'HorizontalAlignment', 'right',...
    'Interpreter', 'tex', 'FontSize', 11, 'FontWeight', 'bold');

% text(real(z*exp(-2i*sigma)), imag(z*exp(-2i*sigma)),...
%     '$f_{\underline X_{n+P} \, |\, \underline X_n}$', 'Interpreter', 'latex', 'FontSize', 12);


% axis square
axis equal
axis([min(x), max(x), min(y), max(y)]);
axis off

fig.Position(3:4) = [350 250];

set(fig,'renderer','Painters');
saveas(fig, 'proba', 'epsc');

