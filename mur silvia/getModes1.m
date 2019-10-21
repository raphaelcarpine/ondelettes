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
Q = 50;
MaxRidges = 10;
MaxParallelRidges = 2;
fmin = 30;
fmax = 40;
NbFreq = 300;

WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges);



%% test

ridges = {};
for k = 1:9
    ridges{end+1} = RidgeExtract(t, X(k,:), Q, fmin, fmax, NbFreq,...
        'NbMaxParallelRidges', 2, 'NbMaxRidges', 10);
end



[t, freqs, shapes, errors, ridgesNumber] = getModes(ridges);


mode = 1;

% figure;
% plot(t{mode}, abs(shapes{mode}));
% ylabel('|T|');
% figure;
% plot(t{mode}, angle(shapes{mode}));
% ylabel('arg(T)');
figure;
plot(t{mode}, real(shapes{mode}));
ylabel('Re(T)');
figure;
plot(t{mode}, imag(shapes{mode}));
ylabel('Im(T)');

figure;
plot(t{mode}, freqs{mode});




shapes0 = nan(1, 9);
for k = 1:9
    shapet = shapes{mode}(k,:);
    shapes0(k) = mean(shapet(~isnan(shapet)));
    if isnan(shapes0(k))
        shapes0(k) = 0;
    end
end

plotModShape(real(shapes0));



%% test complexes

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











