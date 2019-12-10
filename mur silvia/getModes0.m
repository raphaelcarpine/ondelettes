%etape et transient
P = 7;
transient = 0;

t0 = 0;
tf = inf;


[t, X] = getData(P, transient);

X = X(:, t>=t0 & t<tf);
t = t(t>=t0 & t<tf);


X1 = exp(2i*pi*33.9*t) .* exp(- 0.0051*2*pi*33.9*t);

X1 = [0.5587 - 0.0295i
  -0.0026 - 0.0307i
  -0.5400 - 0.0327i
   0.3721 - 0.0224i
  -0.0280 - 0.0249i
  -0.4466 - 0.0248i
   0.2481 - 0.0184i
  -0.0027 - 0.0007i
  -0.0309 + 0.0042i] * X1;

X2 = exp(2i*pi*36.8*t + 1) .* exp(- 0.0065*2*pi*36.8*t);

X2 = [0.6162 - 0.0017i
   0.0447 - 0.0037i
  -0.4815 - 0.0105i
   0.3967 - 0.0064i
   0.0045 + 0.0018i
  -0.4021 + 0.0015i
   0.2585 - 0.0022i
  -0.0002 + 0.0003i
  -0.0287 + 0.0059i] * X2;

%X = real(X1 + X2);



% fig = figure;
% ax = axes(fig);
% hold(ax, 'on');
% plts = nan(1, 9);
% for i = 1:9
%     plts(i) = plot(t, X(i,:), 'Parent', ax);
% end

sensors = 1:9;



fig = figure;
ax = axes(fig);
plts = plot(t, X(sensors,:), 'Parent', ax);
plts = transpose(plts);


%ondelette
Q = 30;
MaxRidges = 1;
MaxParallelRidges = 1;
fmin = 26;
fmax = 30;
NbFreq = 300;
ct = 3;
cf = 5;

WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges, 'CtEdgeEffects', ct);


%% choix de Q

Df = 2.5;
f = 28.5;
Dt = 1.7;

T = t(end) - t(1);

[Qmin, Qmax, Qz] = getBoundsQ(f, Df, Dt, T, ct, cf);

disp(['Qmin = ', num2str(Qmin)]);
disp(['Qmax = ', num2str(Qmax)]);
disp(['Qz = ', num2str(Qz)]);

Q = (Qmin + min(Qmax, Qz)) / 2;
disp(['Q : ', num2str(Q)]);
disp('');


%% test

[time, freq, shape, amplitude] = getModesSingleRidge(t, X, Q, fmin, fmax, NbFreq,...
    'NbMaxRidges', MaxRidges, 'NbMaxParallelRidges', MaxParallelRidges,...
    'ctLeft', ct, 'ctRight', ct);


%%

meanFreq = cell(size(freq));
meanShape = cell(size(shape));
zeta = cell(size(amplitude));
for kr = 1:length(time)
    meanFreq{kr} = mean(freq{kr});
    meanShape{kr} = mean(shape{kr}, 2);
    fig = figure;
    plot(time{kr}, abs(amplitude{kr}));
    Alambda = RegressionMenu('Equation', 'a*exp(-lambda*x)', 'Param', 'a lambda',...
        'Param0', '1 1');
    delete(fig);
    zeta{kr} = Alambda(2) / (2*pi*meanFreq{kr});
end




%%

for mode = 1:length(time)
    
    figure;
    plot(time{mode}, angle(shape{mode}*exp(-1i*pi/2)) + pi/2);
    hold on
    plot(time{mode}, zeros(size(time{mode})), 'black--');
    plot(time{mode}, pi*ones(size(time{mode})), 'black--');
    xlabel('t');
    ylabel('arg(T)');
    figure;
    plot(abs(amplitude{mode}), real(shape{mode}));
    xlabel('|A|');
    ylabel('Re(T)');
    figure;
    plot(abs(amplitude{mode}), imag(shape{mode}));
    xlabel('|A|');
    ylabel('Im(T)');
    figure;
    plot(time{mode}, abs(amplitude{mode}));
    ylabel('|A|');
    
    % figure;
    % plot(t{mode}, real(amplitudes{mode}));
    % ylabel('Re(A)');
    
    figure;
    plot(time{mode}, freq{mode});
    ylabel('f');
    
    
    
    
    
    shape0 = meanShape{mode};
    
    plotModShape(real(shape0), ['freq : ', num2str(meanFreq{mode}),...
        ' ; amort. : ', num2str(zeta{mode})]);
    
    
    
    
    figure;
    for k=1:9
        polarplot([0, shape0(k)], '-o');
        hold on
    end
    
    figure('Position', [400 300 450 450]);
    hold on
    circle = exp(1i*linspace(0, 2*pi, 30));
    circle2 = 1/2 * exp(1i*linspace(0, 2*pi, 30));
    l0 = linspace(-1, 1, 2);
    l30 = exp(1i*pi/6) * l0;
    l60 = exp(1i*pi/3) * l0;
    l90 = exp(1i*pi/2) * l0;
    l120 = exp(2i*pi/3) * l0;
    l150 = exp(5i*pi/6) * l0;
    for k=1:9
        p0 = 2*mod(k-1, 3) + 2*1i*(3-fix((k-1)/3));
        p1 = p0 + shape0(k);
        grey = 0.9;
        plot(real(p0 + circle2), imag(p0 + circle2), 'color', grey+[0 0 0]);
        plot(real(p0 + l0), imag(p0 + l0), 'color', grey+[0 0 0]);
        plot(real(p0 + l30), imag(p0 + l30), 'color', grey+[0 0 0]);
        plot(real(p0 + l60), imag(p0 + l60), 'color', grey+[0 0 0]);
        plot(real(p0 + l90), imag(p0 + l90), 'color', grey+[0 0 0]);
        plot(real(p0 + l120), imag(p0 + l120), 'color', grey+[0 0 0]);
        plot(real(p0 + l150), imag(p0 + l150), 'color', grey+[0 0 0]);
        plot(real(p0 + circle), imag(p0 + circle), 'black');
        plot(real([p0 p1]), imag([p0, p1]), '-o');
        text(real(p0), imag(p0)+0.5, ['ch', num2str(k)]);
    end
    
    pbaspect(gca, [1 1 1]);
    axis off


end


%% enregistrement
% 
% enregistrement = input('enregistrer ?');
% 
% if ~enregistrement
%     return
% end
% 
% title = ['P', num2str(P), 'T', num2str(transient), 





