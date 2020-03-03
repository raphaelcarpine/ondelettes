%etape et transient
P = 0;
transient = 0;

t0 = 1450.2;
tf = 1452;


[t, X] = getData(P, transient);

X = X(:, t>=t0 & t<tf);
t = t(t>=t0 & t<tf);
t = t-t0;



%% affichage optionnel
if false
    
    % fig = figure;
    % ax = axes(fig);
    % hold(ax, 'on');
    % plts = nan(1, 9);
    % for i = 1:9
    %     plts(i) = plot(t, X(i,:), 'Parent', ax);
    % end
    
    sensor = 1;
    
    
    
    fig = figure;
    ax = axes(fig);
    plts = plot(t, X(sensor,:), 'Parent', ax);
    plts = transpose(plts);
    
    
    %ondelette
    Q = 5;
    MaxRidges = 1;
    MaxParallelRidges = 1;
    fmin = 30;
    fmax = 40;
    NbFreq = 100;
    
    WaveletMenu('WaveletPlot', plts, 'fmin', fmin, 'fmax', fmax,...
        'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges);
    
    return
    
end


%% regressions

% wavelet

Q = 5;
% MaxRidges = 1;
% MaxParallelRidges = 1;
fmin = 30;
fmax = 40;
NbFreq = 100;


ridges = cell(1, 9);
for sensor = 1:9
    ridges{sensor} = RidgeExtract(t, X(sensor,:), Q, fmin, fmax, NbFreq,...
        'NbMaxParallelRidges', 1, 'NbMaxRidges', 1);
end


% regressions

paramsReg = nan(9, 8);

paramValue01 = [0.02, 0.01, 1, 1, 3*2*pi, 0]; % a1, a2, l1, l2, dw, delta
paramValue02 = [34*2*pi, 37*2*pi, 0]; % w1, w2, d
paramValue03 = [0, 0]; % delta1, delta2

for sensor = 1:9
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
    
    
    %%%%%%%% reg enveloppe %%%%%%%%
    
    fig = figure;
    plot(tr, abs(Tr));
    xlabel('t');
    ylabel('A_{12}');
    
    
    eq1 = 'sqrt(a1^2*exp(-2*l1*x) + a2^2*exp(-2*l2*x) + 2*a1*a2*exp(-(l1+l2)*x)*cos(dw*x+d))';
    param1 = 'a1 a2 l1 l2 dw d';
    
    paramValue01 = RegressionMenu('Equation', eq1, 'Param', param1, 'Param0', paramValue01);
    
    delete(fig);
    
    a1 = paramValue01(1);
    a2 = paramValue01(2);
    l1 = paramValue01(3);
    l2 = paramValue01(4);
    
    
    %%%%%%%% reg freq %%%%%%%%
    
    fig = figure;
    plot(tr, Fr);
    xlabel('t');
    ylabel('\omega_1 + \dot\psi');
    
    A1 = num2str(a1, 10);
    A2 = num2str(a2, 10);
    A = num2str(a2/a1, 10);
    Dl = num2str(l2-l1, 10);
    
    eq2 = ['w1/(2*pi) + 1/(2*pi)*(-', Dl, '*sin((w2-w1)*x+d) + (w2-w1)*(', A, '*exp(-', Dl, '*x)+cos((w2-w1)*x+d)))',...
        ' / (1/', A, '*exp(', Dl, '*x) + ', A, '*exp(-', Dl, '*x) + 2*cos((w2-w1)*x+d))'];
    param2 = 'w1 w2 d';
    
    
    paramValue02 = RegressionMenu('Equation', eq2, 'Param', param2, 'Param0', paramValue02);
    
    delete(fig);
    
    w1 = paramValue02(1);
    w2 = paramValue02(2);
    d = paramValue02(3);
    
    
    %%%%%%%% reg phase %%%%%%%%
    
    fig = figure;
    plot(tr, angle(Tr.*exp(-1i*w1*tr)));
    xlabel('t');
    ylabel('\delta_1 + \psi');
    
    Dw = num2str(w2-w1, 10);
    
    eq1 = ['d1 + angle(1 + ', A, '*exp(-', Dl, '*x+1i*', Dw, '*x+1i*(d2-d1)))'];
    param1 = 'd1  d2';
    
    paramValue03 = RegressionMenu('Equation', eq1, 'Param', param1, 'Param0', paramValue03);
    
    delete(fig);
    
    
    
    %%%%%%%% reg complexes %%%%%%%%
    
    F = @(p) p(1)*exp(-p(3)*tr+1i*p(5)*tr+1i*p(7)) + p(2)*exp(-p(4)*tr+1i*p(6)*tr+1i*p(8));
    S = @(p) abs(Tr - F(p));
    
%     a0 = 0.01*[2, 0.7];
%     l0 = [0.93, 1.0146];
%     w0 = 2*pi*[33.923, 36.69];
%     d0 = 1.3*[1 0] + 0.05*2*pi*[1 1];
%     p0 = [a0, l0, w0, d0];
%     
%     figure;
%     hold on
%     plot3(tr, imag(Tr), real(Tr), 'b', 'LineWidth', 1);
%     plot3(tr, imag(F(p0)), real(F(p0)), 'r', 'LineWidth', 1);
%     xlabel('t');
%     ylabel('Im');
%     zlabel('Re');
%     view(15, 20);
    
    paramValue04 = [paramValue01(1:4), paramValue02(1:2), paramValue03];
    
    
    % détermination de MaxFunctionEvaluations, pour ne pas calculer plus
    % de 30s
    tic
    for k = 1:100
        S(paramValue04);
    end
    T100 = toc;
    
    maxFunctionEval = round(30*100/T100);
    
    
    optionsReg = optimoptions(@lsqnonlin, 'MaxIterations', inf,...
        'StepTolerance', 1e-10, 'MaxFunctionEvaluations', maxFunctionEval, 'FunctionTolerance', 0);
    
    lb = ones(size(paramValue04))*(-inf);
    ub = ones(size(paramValue04))*inf;
    
    try
        p1 = lsqnonlin(S, paramValue04, lb, ub, optionsReg);
    catch
        p1 = zeros(1, 8);
    end
        
    disp(['a1 = ', num2str(p1(1), 6)]);
    disp(['lambda1 = ', num2str(p1(3), 6)]);
    disp(['w1 = ', num2str(p1(5), 6), ' ; f1 = ', num2str(p1(5)/2/pi, 6)]);
    disp(['delta1 = ', num2str(p1(7), 6)]);
    disp(['a2 = ', num2str(p1(2), 6)]);
    disp(['lambda2 = ', num2str(p1(4), 6)]);
    disp(['w2 = ', num2str(p1(6), 6), ' ; f2 = ', num2str(p1(6)/2/pi, 6)]);
    disp(['delta2 = ', num2str(p1(8), 6)]);
    
    fig = figure;
    hold on
    plot3(tr, imag(Tr), real(Tr), 'b', 'LineWidth', 1);
    plot3(tr, imag(F(p1)), real(F(p1)), 'r', 'LineWidth', 1);
    xlabel('t');
    ylabel('Im');
    zlabel('Re');
    view(15, 20);
    
    paramsReg(sensor, :) = p1;
    
    waitfor(fig);
end


% mean

w1 = mean(paramsReg(:,5));
w2 = mean(paramsReg(:,6));
l1 = mean(paramsReg(:,3));
l2 = mean(paramsReg(:,4));

shapes = {transpose(paramsReg(:,1).*exp(1i*paramsReg(:,7))), transpose(paramsReg(:,2).*exp(1i*paramsReg(:,8)))};

for ks = 1:2
    shapes{ks} = shapes{ks} ./ shapes{ks}(1);
end




%% test

mode = 2;

shapes0 = shapes{mode};

plotModShape(real(shapes0));



%% test complexes

figure;
for k=1:9
    polarplot([0, shapes0(k)], '-o');
    hold on
end

fig = figure;
ax = axes(fig);
hold(ax, 'on');
circle = exp(1i*linspace(0, 2*pi, 30));
plot([-1, 5], [0, 0], '--', 'Color', [0, 0, 0, 0.5]);
plot([-1, 5], [2, 2], '--', 'Color', [0, 0, 0, 0.5]);
plot([-1, 5], [4, 4], '--', 'Color', [0, 0, 0, 0.5]);
plot([0, 0], [-1, 5], '--', 'Color', [0, 0, 0, 0.5]);
plot([2, 2], [-1, 5], '--', 'Color', [0, 0, 0, 0.5]);
plot([4, 4], [-1, 5], '--', 'Color', [0, 0, 0, 0.5]);
for kx = 0:2
    for ky = 0:2
        p0 = 2*kx + 2*1i*(2-ky);
        p1 = p0 + shapes0(3*ky+kx+1);
        plot(real(p0 + circle), imag(p0 + circle), 'black');
        colorI = get(gca,'ColorOrderIndex');
        plot(real([p0 p1]), imag([p0, p1]), 'r', 'LineWidth', 2);
        set(gca,'ColorOrderIndex', colorI);
        plot(real(p1), imag(p1) , 'ro', 'LineWidth', 2);
    end
end

axis equal;
ax.Position = ax.OuterPosition;
% fig.Position = 40*ones(1, 4) + ax.Position;

axis off

pbaspect(gca, [1 1 1]);










