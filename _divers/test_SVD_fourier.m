phi = [ 1,  1,  1;
        2, -1, -1;
        1, -1,  0];
% phi = [ 1,  1, 0;
%         1, -1, 0];
f = [1, 1.02, 15];
z = [0.03, 0.03, 0.02];
d = [0, 0, 0];

phi(:, 3) = 0*phi(:, 3);

%%

L = 2i*pi*f.*sqrt(1-z.^2) - 2*pi*z.*f;

t = linspace(0, 1000, 100000);
x = zeros(size(t));
for k = 1:length(f)
    x = x + real(phi(:, k)*exp(L(k)*t+d(k)));
end

% % test idee
% x = [x; x; x; x];

figure;
plt = plot(t, x);

%%

fmin = 0;
fmax = 20;
lag = 5/(2*pi*f(1)*z(1));

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'RemoveMean', true,...
    'AutocorrelationMode', true, 'AutocorrelationSVDMode', true, 'AutocorrelationFourierSVDMode', true,...
    'AutocorrelationMaxLag', lag, 'AutocorrelationNsvd', 2, 'FourierScale', 'log');
