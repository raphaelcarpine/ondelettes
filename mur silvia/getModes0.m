% close all

%etape et transient
P = 0;
transient = 0;

singleRidgeMode = false;

t0 = 0;
tf = 80;
t0 = 1445;
tf = 1460;
t0 = -inf;
tf = inf;


[t, X] = getData(P, transient);

X = X(:, t>=t0 & t<tf);
t = t(t>=t0 & t<tf);


% rééchantillonnage
indices = 1:3:length(t);
X = X(:, indices);
t = t(indices);


% T0 = 100;
% X1 = exp(2i*pi*33.9*(t-T0)) .* exp(- 0.0051*2*pi*33.9*(t-T0));
% X1(t < T0) = 0;
% shape1 = [0.5587 - 0.0295i;  -0.0026 - 0.0307i;  -0.5400 - 0.0327i;   0.3721 - 0.0224i;  -0.0280 - 0.0249i;  -0.4466 - 0.0248i;   0.2481 - 0.0184i;  -0.0027 - 0.0007i;  -0.0309 + 0.0042i];
% X2 = exp(2i*pi*36.8*(t-T0) + 1) .* exp(- 0.0065*2*pi*36.8*(t-T0));
% X2(t < T0) = 0;
% shape2 = [0.6162 - 0.0017i;   0.0447 - 0.0037i;  -0.4815 - 0.0105i;   0.3967 - 0.0064i;   0.0045 + 0.0018i;  -0.4021 + 0.0015i;   0.2585 - 0.0022i;  -0.0002 + 0.0003i;  -0.0287 + 0.0059i];
% X = real([shape1*X1 + shape2*X2]);


sensors = [1];

%%% test %%%
if false
    for k = 1:9
        x = X(k,:);
        Rx = xcorr(x, 'biased') / var(x);
        Rx = Rx(ceil(length(Rx)/2):end);
        X(k, :) = Rx;
    end
end
%%% test %%%

%%% test %%%
if false
    inputSensor = 1;
    outputSensor = 2;
    
    sensors = 1;
    
    X(sensors, :) = ifft(fft(X(outputSensor, :)) ./ fft(X(inputSensor, :)));
end
%%% test %%%



fig = figure;
ax = axes(fig);
plts = plot(t, X(sensors,:), 'Parent', ax);
xlabel(ax, 'time');
ylabel(ax, 'acceleration');
plts = transpose(plts);


%ondelette
Q = 30;
MaxRidges = 1;
MaxParallelRidges = 1;
fmin = 4;
fmax = 40;
NbFreq = 300;

ct = 3;
cf = 5;

WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax, 'NbFreq', NbFreq,...
    'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges...
    , 'CtEdgeEffects', ct, 'MultiSignalMode', true);


return

%%
try
    close(fig);
catch
    delete(fig);
end



%% choix de Q



Df = 3;
f = 11;
Dt = 1.5;

T = t(end) - t(1);

[Qmin, Qmax, Qz] = getBoundsQ(f, Df, Dt, T, ct, cf);

disp(['Qmin = ', num2str(Qmin)]);
disp(['Qmax = ', num2str(Qmax)]);
disp(['Qz = ', num2str(Qz)]);

Q = (Qmin + min(Qmax, Qz)) / 2;
disp(['Q : ', num2str(Q)]);
disp('');



%% extraction des modes

if singleRidgeMode
    [time, freq, shape, amplitude] = getModesSingleRidge(t, X, Q, fmin, fmax, NbFreq,...
        'NbMaxRidges', MaxRidges, 'NbMaxParallelRidges', MaxParallelRidges,...
        'ctLeft', ct, 'ctRight', ct);
else
    ridges = {};
    for k = 1:9
        ridges{end+1} = RidgeExtract(t, X(k,:), Q, fmin, fmax, NbFreq,...
            'NbMaxParallelRidges', MaxRidges, 'NbMaxRidges', MaxParallelRidges);
    end
    
    [time, freq, freqs, shape, amplitude, errors, ridgesNumber] = getModes(ridges, 1);
end


%%

meanFreq = cell(size(freq));
meanShape = cell(size(shape));
zeta = cell(size(amplitude));
for kr = 1:length(time)
    meanFreq{kr} = mean(freq{kr});
    meanShape{kr} = nan(1, 9);
    for kddl = 1:9
        shapeKddl = shape{kr}(kddl, :);
        meanShape{kr}(kddl) = mean (shapeKddl (~isnan(shapeKddl)));
    end
    fig = figure;
    plot(time{kr}, abs(amplitude{kr}));
    Alambda = RegressionMenu('Equation', 'a*exp(-lambda*x)', 'Param', 'a lambda',...
        'Param0', '1 1');
    delete(fig);
    zeta{kr} = Alambda(2) / (2*pi*meanFreq{kr});
end




%%

for mode = 1:length(time)
    
    figs =  figure;
    plot(time{mode}, (angle(shape{mode}*exp(-1i*pi/2)) + pi/2) * 180/pi);
    hold on
    plot(time{mode}, zeros(size(time{mode})), 'black--');
    plot(time{mode}, 180*ones(size(time{mode})), 'black--');
    ylim([-90, 270]);
    xlabel('t');
    ylabel('arg(T)');
    
    figs(end+1) = figure;
    plot(time{mode}, real(shape{mode}));
    xlabel('t');
    ylabel('Re(T)');
    figs(end+1) = figure;
    plot(time{mode}, imag(shape{mode}));
    xlabel('t');
    ylabel('Im(T)');
    
    figs(end+1) = figure;
    hold on
    plot(time{mode}, abs(amplitude{mode}));
    xlabel('t');
    ylabel('|A|');
    
    figs(end+1) = figure;
    plot(time{mode}, freq{mode});
    xlabel('t');
    ylabel('f');
    
    
    
    
    
    shape0 = meanShape{mode};
    
    plotModShape(real(shape0), ['freq : ', num2str(meanFreq{mode}),...
        ' ; amort. : ', num2str(10*zeta{mode}), '%']);
    
    
    figure;
    for k=1:9
        polarplot([0, shape0(k)], '-o');
        hold on
    end
    
    plotComplexModShape(shape0, ['freq : ', num2str(meanFreq{mode}),...
        ' ; amort. : ', num2str(100*zeta{mode}), '%']);
    
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





