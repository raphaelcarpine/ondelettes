clear all
close all

Transients = 2:3;

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

%% sauvegarde

% nombre de modes pour P0, P6 et P7
nbModes = [3, 4, 3];

% variable de sauvegarde
ModesWvlt2 = struct([]);

%% etape et transient

Transient  = struct([]);

% transient 1
Transient(1).P = 0; % bof en fait
Transient(1).t0 = 1267.35;
Transient(1).tf = 1269.95;
Transient(1).f1 = 33.93;
Transient(1).f2 = 36.76;
Transient(1).z1 = 0.0029;
Transient(1).z2 = 0.0056;
Transient(1).f3 = 8.35;
Transient(1).noiseCoeff = 0;

% transient 2
Transient(2).P = 0;
Transient(2).t0 = 1450.2; %1450
Transient(2).tf = 1455; %1452
Transient(2).f1 = 33.93;
Transient(2).f2 = 36.76;
Transient(2).z1 = 0.0029;
Transient(2).z2 = 0.0056;
Transient(2).f3 = 8.35;
Transient(2).noiseCoeff = 5;

% transient 3
Transient(3).P = 7;

% Transient(3).t0 = 36.81; % non
% Transient(3).tf = 39;

% Transient(3).t0 = 61.7; % non
% Transient(3).tf = 63.2;

Transient(3).t0 = 305.9; % oui sans threshold, non avec
Transient(3).tf = 306.8;
Transient(3).tf = 308;

% Transient(3).t0 = 458.1; % non
% Transient(3).tf = 459.8;

Transient(3).f1 = 28.23;
Transient(3).f2 = 34.15;
Transient(3).z1 = 0.0066;
Transient(3).z2 = 0.0070;
Transient(3).f3 = 11.02;
Transient(3).noiseCoeff = 2.5;

for kTransient = Transients
    
    %% initialisation
    P = Transient(kTransient).P;
    t0 = Transient(kTransient).t0;
    tf = Transient(kTransient).tf;
    f1 = Transient(kTransient).f1;
    f2 = Transient(kTransient).f2;
    lambda1 = Transient(kTransient).z1 * 2*pi*f1;
    lambda2 = Transient(kTransient).z2 * 2*pi*f2;
    f3 = Transient(kTransient).f3;
    noiseCoeff = Transient(kTransient).noiseCoeff;
    
    [t, X] = getData(P, 0);
    
    %% calcul de Q
    [Qmin, Qmax, Qa] = getBoundsQsingleRidge(f1, f2, lambda1, lambda2, f3, tf-t0, ct, cf);
    if Qmin > Qa
        warning('Qmin > Qa');
    end
    if Qmin > Qmax
        warning('Qmin > Qmax');
    end
    
    Q = (Qmin + Qa)/2;
    
    disp(['Qmin = ', num2str(Qmin), ' ; Qa = ', num2str(Qa), ' ; Qmax = ', num2str(Qmax)]);
    disp(['Q = ', num2str(Q)]);
    
    
    %% affichage optionnel
    if showTest
        sensor = 1:9;
        
        fig = figure;
        ax = axes(fig);
        plts = plot(t, X(sensor,:), 'Parent', ax);
        plts = transpose(plts);
        
        
        %ondelette
        %     Q = 5;
        MaxRidges = 1;
        MaxParallelRidges = 1;
        fmin = 20;
        fmax = 45;
        NbFreq = 300;
        
        WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax,...
            'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges,...
            'XLim', [t0, tf]);
        
        str = input('continue ? ', 's');
        
        if isequal(str, '0') || isequal(str, 'n') || isequal(str, 'no') || isequal(str, 'non')
            return
        end
    end
    
    X = X(:, t>=t0 & t<tf);
    t = t(t>=t0 & t<tf);
    t = t-t0;
    
    
    
    %% ridge extraction
    
    % wavelet
    fmin = 20;
    fmax = 45;
    NbFreq = 300;
    
    
    ridges = cell(1, 9);
    for sensor = 1:9
        ridges{sensor} = RidgeExtract(t, X(sensor,:), Q, fmin, fmax, NbFreq,...
            'NbMaxParallelRidges', 1, 'NbMaxRidges', 1);
    end
    
    
    %% regressions
    
    paramsReg = nan(9, 4);
    
    if P == 0
        paramValue10 = [0.02, 0.01, 1, 1, 2.9*2*pi, 0]; % a1, a2, l1, l2, dw, delta
        paramValue20 = [33.9*2*pi, 36.8*2*pi, 0]; % w1, w2, d
        paramValue30 = [0, 0]; % delta1, delta2, dw
    elseif P == 7
        paramValue10 = [0.02, 0.01, 1, 1, 6*2*pi, 0]; % a1, a2, l1, l2, dw, delta
        paramValue20 = [34.2*2*pi, 28.3*2*pi, 0]; % w1, w2, d
        paramValue30 = [0, 0]; % delta1, delta2
    end
    
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
    
    for sensor = 1:9
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
        
        [t, X] = getData(P, 0);
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
    
    
    
    
    %% affichage
    
    for mode = 1:2
        plotModShape(real(shapes{mode}), ['P', num2str(P), 'M', num2str(mode+1)]);
        plotComplexModShape(shapes{mode}, ['P', num2str(P), 'M', num2str(mode+1), '_complex']);
    end
    
    %% sauvegarde
    for mode = 1:2
        if mode == 1
            p = p1;
        else
            p = p2;
        end
        kP = find([0 6 7] == P);
        ModesWvlt2(kP, mode+1).freq = imag(p)/2/pi;
        ModesWvlt2(kP, mode+1).damping = -real(p)/imag(p);
        ModesWvlt2(kP, mode+1).shape = shapes{mode};
        
        ModesWvlt2(kP, mode+1).choixQ = ['Q = ', num2str(Q), ' (Qmin = ', num2str(Qmin),...
            ' ; Qa = ', num2str(Qa), ' ; Qmax = ', num2str(Qmax), ')'];
        ModesWvlt2(kP, mode+1).choixNoise = ['threshold = ', num2str(noiseCoeff), '*', char(noiseMean)];
    end
    
    
    
end


%% sauvegarde

if saveFiles
    directory = 'mur silvia\modesOndelette2\';
    save([directory, 'ModesWvlt2.mat'], 'ModesWvlt2');
end




