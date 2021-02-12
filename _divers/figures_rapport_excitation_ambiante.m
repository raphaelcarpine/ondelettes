%% complexes Xn Xn+k

z0 = 0.6*exp(0.9i);
zeta = 0.01;
r0 = exp(-2*pi*zeta/sqrt(1-zeta^2));
z = r0*z0;
sigma = abs(z0) * 2*sqrt(pi*zeta/sqrt(1-zeta^2));

x = linspace(-0.5, 1, 300);
y = linspace(-0.5, 1, 300);

[X, Y] = meshgrid(x, y);


P = -exp(-((X-real(z)).^2 + (Y-imag(z)).^2)/(2*sigma^2));

figure;
p = pcolor(X, Y, P);
p.FaceColor = 'interp';
p.EdgeColor = 'none';
colormap gray %bone

yline(0, '-', 'Re');
xline(0, '-', 'Im', 'LabelHorizontalAlignment', 'left', 'LabelOrientation', 'horizontal');

hold on
plot(real(z0), imag(z0), 'r+', 'LineWidth', 1.5, 'MarkerSize', 9);
plot(real(z), imag(z), 'g.', 'LineWidth', 1, 'MarkerSize', 7);

text(real(z0), imag(z0), {'  $\underline X_n$', ''}, 'Color', 'r', 'Interpreter', 'latex', 'FontSize', 12);
text(real(z), imag(z), {'', '$r_0\underline X_n$'}, 'Color', 'g',...
    'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'FontSize', 12);

text(real(z*exp(-2i*sigma)), imag(z*exp(-2i*sigma)),...
    '$f_{\underline X_{n+K} \, |\, \underline X_n}$', 'Interpreter', 'latex', 'FontSize', 12);


axis([-0.5 1 -0.5 1]);
axis square
axis off