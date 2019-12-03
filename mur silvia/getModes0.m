%etape et transient
P = 6;
transient = 1;

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

sensors = 1:9;



fig = figure;
ax = axes(fig);
plts = plot(t, X(sensors,:), 'Parent', ax);
plts = transpose(plts);


%ondelette
Q = 30;
MaxRidges = 2;
MaxParallelRidges = 2;
fmin = 7;
fmax = 13;
NbFreq = 300;

WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges);



%% test

[t, freq, shapes, amplitudes] = getModesSingleRidge(t, X, Q, fmin, fmax, NbFreq,...
    'NbMaxRidges', MaxRidges, 'NbMaxParallelRidges', MaxParallelRidges);

%%

for mode = 1:length(t)
    
    figure;
    plot(t{mode}, angle(shapes{mode}*exp(-1i*pi/2)) + pi/2);
    hold on
    plot(t{mode}, zeros(size(t{mode})), 'black--');
    plot(t{mode}, pi*ones(size(t{mode})), 'black--');
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
    plot(t{mode}, abs(amplitudes{mode}));
    ylabel('|A|');
    
    % figure;
    % plot(t{mode}, real(amplitudes{mode}));
    % ylabel('Re(A)');
    
    figure;
    plot(t{mode}, freq{mode});
    ylabel('f');
    
    
    
    
    
    shapes0 = nan(1, 9);
    for k = 1:9
        shapet = shapes{mode}(k,:);
        shapes0(k) = mean(shapet(~isnan(shapet)));
        if isnan(shapes0(k))
            shapes0(k) = 0;
        end
    end
    
    plotModShape(real(shapes0));
    
    
    
    
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








