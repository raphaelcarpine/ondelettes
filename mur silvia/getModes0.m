%etape et transient
P = 0;
transient = 3;

t0 = 0;
tf = inf;


[t, X] = getData(P, transient);

X = X(:, t>=t0 & t<tf);
t = t(t>=t0 & t<tf);





% fig = figure;
% ax = axes(fig);
% hold(ax, 'on');
% plts = nan(1, 9);
% for i = 1:9
%     plts(i) = plot(t, X(i,:), 'Parent', ax);
% end

sensors = 1;



fig = figure;
ax = axes(fig);
plts = plot(t, X(sensors,:), 'Parent', ax);
plts = transpose(plts);


%ondelette
Q =30;
MaxRidges = 2;
MaxParallelRidges = 2;
fmin = 32;
fmax = 46;
NbFreq = 300;

WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges);


%% choix de Q

Df = 3;
f = 35.5;
Dt = 1.7;

T = t(end) - t(1);

[Qmin, Qmax, Qz] = getBoundsQ(f, Df, Dt, T);

disp(['Qmin = ', num2str(Qmin)]);
disp(['Qmax = ', num2str(Qmax)]);
disp(['Qz = ', num2str(Qz)]);
disp('');

Q = (Qmin + min(Qmax, Qz)) / 2;


%% test

[time, freq, shape, amplitude] = getModesSingleRidge(t, X, Q, fmin, fmax, NbFreq,...
    'NbMaxRidges', MaxRidges, 'NbMaxParallelRidges', MaxParallelRidges);


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
    plot(time{mode}, angle(shapes{mode}*exp(-1i*pi/2)) + pi/2);
    hold on
    plot(time{mode}, zeros(size(time{mode})), 'black--');
    plot(time{mode}, pi*ones(size(time{mode})), 'black--');
    xlabel('t');
    ylabel('arg(T)');
    figure;
    plot(abs(amplitudes{mode}), real(shapes{mode}));
    xlabel('|A|');
    ylabel('Re(T)');
    figure;
    plot(abs(amplitudes{mode}), imag(shapes{mode}));
    xlabel('|A|');
    ylabel('Im(T)');
    figure;
    plot(time{mode}, abs(amplitudes{mode}));
    ylabel('|A|');
    
    % figure;
    % plot(t{mode}, real(amplitudes{mode}));
    % ylabel('Re(A)');
    
    figure;
    plot(time{mode}, freq{mode});
    ylabel('f');
    
    
    
    
    
    shapes0 = meanShape{mode};
    
    plotModShape(real(shapes0), ['freq : ', num2str(freq{mode}),...
        ' ; amort. : ', num2str(zeta{mode})]);
    
    
    
    
    figure;
    for k=1:9
        polarplot([0, shapes0(k)], '-o');
        hold on
    end
    
    figure;
    hold on
    circle = exp(1i*linspace(0, 2*pi, 30));
    for k=1:9
        p0 = 2*mod(k-1, 3) + 2*1i*(3-fix((k-1)/3));
        p1 = p0 + shapes0(k);
        plot(real(p0 + circle), imag(p0 + circle), 'black');
        plot(real([p0 p1]), imag([p0, p1]), '-o');
    end
    
    pbaspect(gca, [1 1 1]);


end








