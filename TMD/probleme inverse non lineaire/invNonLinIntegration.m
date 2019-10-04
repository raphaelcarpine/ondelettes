nIalpha = 1000;


mu = 0.01;
w0 = 2*pi;
w1 = 2*pi/(1+mu);
epsilon = 0.1*2*w1*sqrt(3*mu/8/(1+mu));
alpha = 1.5;




T = 200;
dt = 1e-2;

xi = [0; 0];
vi = [1; 0];




%% integration temporelle


nT = 100;

Xi = xi;
Vi = vi;

Dtemp = @(t, Y) [
    Y(3);
    Y(4);
    -w0^2*Y(1) + mu*w1^2*(Y(2)-Y(1)) + mu*epsilon*sign(Y(4)-Y(3))*abs(Y(4)-Y(3))^alpha;
    -w1^2*(Y(2)-Y(1))-epsilon*sign(Y(4)-Y(3))*abs(Y(4)-Y(3))^alpha
    ];



options = odeset('RelTol', 1e-10, 'Stats', 'off', 'MaxStep', 1/(w0*nT));

[t0, Y] = ode45(Dtemp, [0 T], [Xi; Vi], options);


xtemp = Y(:,1);
xtemp2 = Y(:,2);

t2 = linspace(t0(1), t0(end), T*w0/(2*pi)*nT);
xtemp = interp1(t0, xtemp, t2);
xtemp2 = interp1(t0, xtemp2, t2);

t0 = t2;



fig = figure;
ax = axes(fig);
hold(ax, 'on');
waveletplot = plot(t0, xtemp, 'Parent', ax);
waveletplot2 = plot(t0, xtemp2, 'Parent', ax);
uistack(waveletplot2, 'bottom');
grid(ax, 'on');
xlabel(ax, '$t$');
ylabel(ax, '$X$');





%% ondelette


Q = 30;
MaxRidges = 2;
MaxParallelRidges = 2;
fmin = 0.9;
fmax = 1.1;
NbFreq = 100;


WaveletMenu('WaveletPlot', waveletplot, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges);

WaveletMenu('WaveletPlot', waveletplot2, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges);



RegressionMenu;

[t01, At1] = PlotExtract;
[t02, At2] = PlotExtract;

[t012, At12] = PlotExtract;
[t022, At22] = PlotExtract;



%% inversion
om1 = (0.94654 + 0.94667)*pi;
om2 = (1.046 + 1.0459)*pi;

%%
t011 = linspace(max(t01(1), t012(1)), min(t01(end), t012(end)), length(t01));
figure;
plot(t011, interp1(t01, At1, t011) ./ interp1(t012, At12, t011));

t021 = linspace(max(t02(1), t022(1)), min(t02(end), t022(end)), length(t02));
figure;
plot(t021, interp1(t02, At2, t021) ./ interp1(t022, At22, t021));

%%
Om21 = om1 / sqrt(1 - 0.08601);
Om22 = om2 / sqrt(1 + 0.11597);
Om2 = (Om21 + Om22)/2;

Om1 = sqrt((om1^2+om2^2)^2-(om1^2-om2^2)^2) / (2*Om2);
mu = (om1^2+om2^2-Om1^2-Om2^2)/Om2^2;


%%



At1 = abs (sqrt((Om2^2-om1^2)^2+mu*Om2^4)/(Om2^2-om1^2)) * At1;
At2 = abs (sqrt((Om2^2-om2^2)^2+mu*Om2^4)/(Om2^2-om2^2)) * At2;

At12 = abs (sqrt((Om2^2-om1^2)^2+mu*Om2^4)/(Om2^2)) * At12;
At22 = abs (sqrt((Om2^2-om2^2)^2+mu*Om2^4)/(Om2^2)) * At22;


%%
figure;
hold on
plot(t01, At1);
plot(t02, At2);
% plot(t012, At12);
% plot(t022, At22);




%% nb de points regression
%nPointsReg = 1000;



%% test

S = @(Param) S0(Param, mu, om1, om2, Om1, Om2, xi, vi, t01, t02, nIalpha, At1, At2);

param1 = linspace(0, 0.15, 30);
param2 = linspace(1, 2, 30);
[param13D, param23D] = meshgrid(param1, param2);

Z = nan(length(param1), length(param2));
for ind1 = 1:length(param1)
    for ind2 = 1:length(param2)
        Z(ind1,ind2) = sum (S([param1(ind1), param2(ind2)]).^2);
    end
end

figure;
surf(param13D, param23D, log(Z));
xlabel('epsilon');
ylabel('alpha');





%% regression




Param0 = [0.1, 1.5];
Param0 = [epsilon, alpha];

maxIter = 1e4;
paramPrecision = 1e-6;

optionsReg = optimoptions(@lsqnonlin, 'MaxIterations', maxIter,...
   'StepTolerance', paramPrecision, 'MaxFunctionEvaluations', inf, 'FunctionTolerance', 0);

Param = lsqnonlin(S, Param0, [], [], optionsReg)


%% courbes regression

[t, A1, A2] = ridgesAprox(mu, om1, om2, Om1, Om2, Param(1), Param(2), xi, vi, t01, t02, nIalpha);

figure;
hold on
plot(t01, At1, 'b');
plot(t02, At2, 'b');
plot(t, A1, 'r--');
plot(t, A2, 'r--');
xlabel('$t$');
ylabel('$A$');




%% fonctions





    function [A1, A2] = ridgesAprox(mu, om1, om2, Om1, Om2, epsilon, alpha, x0, v0, t01, t02, nIalpha)
        % paramètres
        lambda1 = om1^2 / sqrt((om1^2-Om2^2)^2+mu*Om2^4);
        lambda2 = om2^2 / sqrt((om2^2-Om2^2)^2+mu*Om2^4);
        D1 = lambda1;
        D2 = lambda2;
        
        G1 = 1/2*D1*lambda1^alpha*om1^(alpha-1)*pi^(-3/2)*gamma(alpha/2+1)/gamma(alpha/2+3/2);
        G2 = 1/2*D2*lambda2^alpha*om2^(alpha-1)*pi^(-3/2)*gamma(alpha/2+1)/gamma(alpha/2+3/2);
        
        
        matM12 = diag([1 sqrt(mu)]);
        matO = [(Om2^2-om1^2)/sqrt((om1^2-Om2^2)^2+mu*Om2^4), (Om2^2-om2^2)/sqrt((om2^2-Om2^2)^2+mu*Om2^4);...
            sqrt(mu)*Om2^2/sqrt((om1^2-Om2^2)^2+mu*Om2^4), sqrt(mu)*Om2^2/sqrt((om2^2-Om2^2)^2+mu*Om2^4)];
        
        % conditions initiales
        X0 = x0;
        V0 = v0;
        
        Za0 = 1i*transpose(matO)*matM12*X0 + diag([om1 om2])\transpose(matO)*matM12*V0;
        A0 = abs(Za0);
        
        % integration
        D = @(t, A) -epsilon*mu * [G1 * Ialpha(lambda2*om2*A(2)/(lambda1*om1*A(1)), alpha, nIalpha) * A(1)^alpha ;...
            G2 * Ialpha(lambda1*om1*A(1)/(lambda2*om2*A(2)), alpha, nIalpha) * A(2)^alpha];
        
        
        
        options0 = odeset('RelTol', 1e-10, 'Stats', 'off'); %, 'MaxStep', 1/(w0*nT)
        
        [t, A] = ode45(D, [0, max(t01(end), t02(end))], A0, options0);
        
        A = transpose(A);
        
        A1 = interp1(t, A(1,:), t01);
        A2 = interp1(t, A(2,:), t02);
    end

    
    
    
    function I = Ialpha(lambda, alpha, nIalpha)
        if alpha<1 && abs(lambda) == inf
            I = 0;
            return
        end
        theta = linspace(0, 2*pi, nIalpha);
        theta = theta(1:end-1);
        I = (1 + lambda^2 + 2*lambda*cos(theta)).^(alpha/2-1/2) .* (1 + lambda*cos(theta));
        I = sum(I)*(theta(2)-theta(1));

%         I = 2*pi * (1+lambda.^2).^((alpha-1)/2) * (1 + (alpha-1)/2*lambda^2/(1+lambda^2));
        
        I = 2*pi;
    end

    
    
    



    function deltax = S0(Param, mu, om1, om2, Om1, Om2, xi, vi, t01, t02, nIalpha, At1, At2)
        Pepsilon = Param(1);
        Palpha = Param(2);
        
%         Pmu = Param(1);
%         Pepsilon = Param(2);
%         Palpha = Param(3);
%         Px0 = Param(4:5);
%         Pv0 = Param(6:7);
        
%         Om2t = sqrt(1/2/(1+Pmu) * (c1 - sqrt(c1^2-(1+Pmu)*c2)));
%         Om1t = sqrt(c2)/2/Om2t;
        
        [A1, A2] = ridgesAprox(mu, om1, om2, Om1, Om2, Pepsilon, Palpha, xi, vi, t01, t02, nIalpha);
        deltax = [A1 - At1, A2 - At2];
    end


    
    
    
    
    



