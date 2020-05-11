t; % vecteur tps
X;
nbSensors = size(X, 1);

%%

% choix precision
ct = 3;
cf = 5;

% choix bruit
noiseMean = @(x) exp( mean( log(x)));
% noiseMean = @mean;

% affichage
showTest = false;
dashedLines = false;
setBounds = false;
saveFiles = false;

%% transient

f1 = 10;
f2 = 10;
z1 = 0.01;
z2 = 0.01;
f3 = 20;
noiseCoeff = 1;

lambda1 = z1 * 2*pi*f1;
lambda2 = z2 * 2*pi*f2;



%% calcul de Q
[Qmin, Qmax, Qa] = getBoundsQsingleRidge(f1, f2, lambda1, lambda2, f3, T, ct, cf);
if Qmin > Qa
    warning('Qmin > Qa');
end
if Qmin > Qmax
    warning('Qmin > Qmax');
end

Q = (Qmin + Qa)/2;
Q = Qmin;

disp(['Qmin = ', num2str(Qmin), ' ; Qa = ', num2str(Qa), ' ; Qmax = ', num2str(Qmax)]);
disp(['Q = ', num2str(Q)]);


%% affichage optionnel
if showTest
    sensor = 1:nbSensors;
    
    fig = figure;
    ax = axes(fig);
    plts = plot(t, X(sensor,:), 'Parent', ax);
    plts = transpose(plts);
    
    
    %ondelette
    %     Q = 5;
    MaxRidges = 1;
    MaxParallelRidges = 1;
    fmin = 5;
    fmax = 25;
    NbFreq = 300;
    
    WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax,...
        'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges,...
        'XLim', [0, T]);
    
    str = input('continue ? ', 's');
    
    if isequal(str, '0') || isequal(str, 'n') || isequal(str, 'no') || isequal(str, 'non')
        return
    end
end



%% ridge extraction

% wavelet
fmin = 5;
fmax = 15;
NbFreq = 300;


ridges = cell(1, nbSensors);
for sensor = 1:nbSensors
    ridges{sensor} = RidgeExtract(t, X(sensor,:), Q, fmin, fmax, NbFreq,...
        'NbMaxParallelRidges', 1, 'NbMaxRidges', 1);
end


%% regressions

paramsReg = nan(nbSensors, 4);

paramValue10 = [0.0002, 0.0001, 1, 1, 0.01, 0]; % a1, a2, l1, l2, dw, delta
paramValue20 = [10*2*pi, 10.01*2*pi, 0]; % w1, w2, d
paramValue30 = [0, 0]; % delta1, delta2, dw

if setBounds
    paramValue1bounds = [0, 0, 0, 0, 0, -inf;
        inf, inf, 0.02*35*2*pi, 0.02*35*2*pi, 10*2*pi, inf];
    paramValue2bounds = [0, 0, -inf;
        50*2*pi, 50*2*pi, inf];
    paramValue3bounds = [-inf, -inf;
        inf, inf];
else
    paramValue1bounds = nan(2, 0);
    paramValue2bounds = nan(2, 0);
    paramValue3bounds = nan(2, 0);
end

for sensor = 1
    disp(['ch', num2str(sensor)]);
    
    ridge = ridges{sensor};
    
    tr = ridge.time{1};
    Tr = ridge.val{1};
    Fr = ridge.freq{1};
    
    Tr = Tr(~isnan(tr));
    Fr = Fr(~isnan(tr));
    tr = tr(~isnan(tr));
    Fr = Fr(~isnan(Tr));
    tr = tr(~isnan(Tr));
    Tr = Tr(~isnan(Tr));
    tr = tr(~isnan(Fr));
    Tr = Tr(~isnan(Fr));
    Fr = Fr(~isnan(Fr));
    
    
    %% %%%%% noise %%%%%%%%
    
    % threshold noise
    
    WvltF = WvltComp(t, X(sensor, :), (f1+f2)/2, Q);
    thresholdNoise = noiseCoeff * noiseMean( abs( WvltF));
    
    %     figure;
    %     plot(t, abs(WvltF));
    %     hold on
    %     plot(t, thresholdNoise * ones(size(t)));
    %     hold off;
    
    % time threshold
    kTfThreshold = 1;
    while kTfThreshold <= length(Tr) && abs(Tr(kTfThreshold)) > thresholdNoise
        kTfThreshold = kTfThreshold + 1;
    end
    if kTfThreshold > length(Tr)
        warning('threshold non atteint');
    elseif kTfThreshold == 1
        warning('threshold non dépassé');
        kTfThreshold = kTfThreshold + 1;
    end
    kTfThreshold = kTfThreshold - 1;
    
    TfThreshold = tr(kTfThreshold);
    
    % threshold
    timeThreshold = tr(1:kTfThreshold);
    freqThreshold = Fr(1:kTfThreshold);
    amplitudeThreshold = Tr(1:kTfThreshold);
    
    
    %% %%%%% reg enveloppe %%%%%%%%
    
    fig = figure('Name', ['ch', num2str(sensor)]);
    ax = axes(fig);
    hold on
    ax.ColorOrderIndex = 2;
    plot(tr, thresholdNoise * ones(size(tr)), 'LineWidth', 1);
    ax.ColorOrderIndex = 1;
    plot(tr, abs(Tr), ':', 'LineWidth', 1);
    ax.ColorOrderIndex = 1;
    plot(timeThreshold, abs(amplitudeThreshold), 'LineWidth', 1);
    hold off
    xlabel('t');
    ylabel('A_{12}');
    
    
    eq1 = 'sqrt(a1^2*exp(-2*l1*x) + a2^2*exp(-2*l2*x) + 2*a1*a2*exp(-(l1+l2)*x)*cos(dw*x+d))';
    param1 = 'a1 a2 l1 l2 dw d';
    
    paramValue1 = RegressionMenu('Equation', eq1, 'Param', param1, 'Param0', paramValue10,...
        'LowerBoundsParams', paramValue1bounds(1, :), 'UpperBoundsParams', paramValue1bounds(2, :));
    
    delete(fig);
    
    if isempty(paramValue1)
        paramsReg(sensor, :) = [0, 0, nan, nan];
        continue
    end
    
    
    a1 = paramValue1(1);
    a2 = paramValue1(2);
    l1 = paramValue1(3);
    l2 = paramValue1(4);
    
    
    %% %%%%% reg freq %%%%%%%%
    
    fig = figure('Name', ['ch', num2str(sensor)]);
    ax = axes(fig);
    hold on
    if dashedLines
        plot(tr, Fr, ':', 'LineWidth', 1);
        ax.ColorOrderIndex = 1;
    end
    plot(timeThreshold, freqThreshold, 'LineWidth', 1);
    hold off
    xlabel('t');
    ylabel('\omega_1 + \psi''');
    
    A1 = num2str(a1, 10);
    A2 = num2str(a2, 10);
    A = num2str(a2/a1, 10);
    Dl = num2str(l2-l1, 10);
    
    %paramValue20(2) = paramValue20(1) + paramValue10(5);
    
    eq2 = ['w1/(2*pi) + 1/(2*pi)*(-', Dl, '*sin((w2-w1)*x+d) + (w2-w1)*(', A, '*exp(-', Dl, '*x)+cos((w2-w1)*x+d)))',...
        ' / (1/', A, '*exp(', Dl, '*x) + ', A, '*exp(-', Dl, '*x) + 2*cos((w2-w1)*x+d))'];
    param2 = 'w1 w2 d';
    
    
    paramValue2(1, :) = RegressionMenu('Equation', eq2, 'Param', param2, 'Param0', paramValue20,...
        'LowerBoundsParams', paramValue2bounds(1, :), 'UpperBoundsParams', paramValue2bounds(2, :));
    paramValue2(2, :) = RegressionMenu('Equation', eq2, 'Param', param2, 'Param0', paramValue20([2, 1, 3]),...
        'LowerBoundsParams', paramValue2bounds(1, :), 'UpperBoundsParams', paramValue2bounds(2, :));
    
    delete(fig);
    
    if norm(sort(paramValue2(1, 1:2)) - sort(paramValue20(1:2)))...
            < norm(sort(paramValue2(2, 1:2)) - sort(paramValue20(1:2)))
        paramValue2 = paramValue2(1, :);
        disp('1');
    else
        paramValue2 = paramValue2(2, :);
        disp('2');
    end
    
    w1 = paramValue2(1);
    w2 = paramValue2(2);
    d = paramValue2(3);
    
    %paramValue01(5) = w2 - w1;
    
    %% %%%%% reg phase %%%%%%%%
    
    phaseShift = angle(Tr.*exp(-1i*w1*tr));
    phaseShift = phaseContinuity(phaseShift);
    phaseShiftThreshold = angle(amplitudeThreshold.*exp(-1i*w1*timeThreshold));
    phaseShiftThreshold = phaseContinuity(phaseShiftThreshold);
    
    fig = figure('Name', ['ch', num2str(sensor)]);
    ax = axes(fig);
    hold on
    if dashedLines
        plot(tr, phaseShift, ':', 'LineWidth', 1);
        ax.ColorOrderIndex = 1;
    end
    plot(timeThreshold, phaseShiftThreshold, 'LineWidth', 1);
    hold off
    xlabel('t');
    ylabel('\delta_1 + \psi');
    
    Dw = num2str(w2-w1, 10);
    
    eq3 = ['d1 + angle(1 + ', A, '*exp(-', Dl, '*x+1i*', Dw, '*x+1i*(d2-d1)))'];
    eq3 = ['phaseContinuity(', eq3, ')'];
    param3 = 'd1  d2';
    
    paramValue3 = RegressionMenu('Equation', eq3, 'Param', param3, 'Param0', paramValue30,...
        'LowerBoundsParams', paramValue3bounds(1, :), 'UpperBoundsParams', paramValue3bounds(2, :));
    
    delete(fig);
    
    
    
    %% %%%%% reg complexes %%%%%%%%
    
    fig = figure('Name', ['ch', num2str(sensor)]);
    ax = axes(fig);
    hold on
    if dashedLines
        plot3(tr, real(Tr), imag(Tr), ':', 'LineWidth', 1);
        ax.ColorOrderIndex = 1;
    end
    plot3(timeThreshold, real(amplitudeThreshold), imag(amplitudeThreshold), 'LineWidth', 1);
    hold off
    xlabel('t');
    ylabel('Re');
    zlabel('Im');
    if dashedLines
        xlim([tr(1), tr(end)]);
    else
        xlim([timeThreshold(1), timeThreshold(end)]);
    end
    view(15, 20);
    
    eq4 = 'C1*exp(p1*x) + C2*exp(p2*x)';
    param4 = 'C1  C2  p1  p2';
    
    w1 = paramValue2(1);
    w2 = paramValue2(2);
    lambda1 = paramValue1(3);
    lambda2 = paramValue1(4);
    a1 = paramValue1(1);
    a2 = paramValue1(2);
    d1 = paramValue3(1);
    d2 = paramValue3(2);
    
    paramValue40 = [a1*exp(1i*d1), a2*exp(1i*d2), -lambda1+1i*w1, -lambda2+1i*w2];
    
    
    % détermination de MaxFunctionEvaluations, pour ne pas calculer plus
    % de 30s
    F = @(p) p(1)*exp(p(3)*tr) + p(2)*exp(p(4)*tr);
    S = @(p) abs(Tr - F(p));
    tic
    for k = 1:100
        S(paramValue40);
    end
    T100 = toc;
    maxFunctionEval = round(30*100/T100);
    
    
    optionsReg = optimoptions(@lsqnonlin, 'MaxIterations', inf,...
        'StepTolerance', 1e-10, 'MaxFunctionEvaluations', maxFunctionEval, 'FunctionTolerance', 0);
    
    try
        paramValue4 = RegressionMenu('Equation', eq4, 'Param', param4, 'Param0', paramValue40,...
            'OptionsRegression', optionsReg);
    catch
        paramValue4 = zeros(1, 4);
    end
    
    delete(fig);
    
    if imag(paramValue4(3)) > imag(paramValue4(4)) % f1 < f2
        paramValue4 = paramValue4([2, 1, 4, 3]);
    end
    
    %% affichage
    
    disp(['f1 = ', num2str(imag(paramValue4(3))/2/pi, 6)]);
    disp(['f2 = ', num2str(imag(paramValue4(4))/2/pi, 6)]);
    disp(['z1 = ', num2str(-real(paramValue4(3))/imag(paramValue4(3))*100, 6), '%']);
    disp(['z2 = ', num2str(-real(paramValue4(4))/imag(paramValue4(4))*100, 6), '%']);
    disp(['a1 = ', num2str(abs(paramValue4(1)), 6)]);
    disp(['a2 = ', num2str(abs(paramValue4(2)), 6)]);
    disp(['delta1 = ', num2str(angle(paramValue4(1)), 6)]);
    disp(['delta2 = ', num2str(angle(paramValue4(2)), 6)]);
    
    paramsReg(sensor, :) = paramValue4;
end


%% mean

C1 = paramsReg(:,1);
C2 = paramsReg(:,2);

p1 = paramsReg(:,3);
p2 = paramsReg(:,4);
p1 = wmean(p1, abs(C1));
p2 = wmean(p2, abs(C2));

disp(' ')
disp('%%%%% mean %%%%%')
disp(' ')
disp(['f1 = ', num2str(imag(p1)/2/pi, 6)]);
disp(['f2 = ', num2str(imag(p2)/2/pi, 6)]);
disp(['z1 = ', num2str(-real(p1)/imag(p1)*100, 6), '%']);
disp(['z2 = ', num2str(-real(p2)/imag(p2)*100, 6), '%']);

shapes = {C1, C2};

for ks = 1:2
    shapes{ks} = shapes{ks} / sqrt(shapes{ks}.' * shapes{ks});
    shapes{ks} = shapes{ks} * sign(real(shapes{ks}(1)));
end







