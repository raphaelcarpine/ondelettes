function invNonLin

nIalpha = 1000;


mu = 0.01;
w0 = 2*pi;
w1 = 2*pi/(1+mu);
epsilon = 0.05*2*w1*sqrt(3*mu/8/(1+mu));
alpha = 1.5;




T = 100;
dt = 1e-2;

xi = [0; 0];
vi = [1; 0];

%% dérivation


    function DA = derivAprox(mu, w0, w1, epsilon, alpha, A)
        % paramètres
        Om1 = w0;
        Om2 = w1;
        om1 = sqrt(1/2 * (Om1^2+(1+mu)*Om2^2 - sqrt((Om1^2+(1+mu)*Om2^2)^2-4*Om1^2*Om2^2)));
        om2 = sqrt(1/2 * (Om1^2+(1+mu)*Om2^2 + sqrt((Om1^2+(1+mu)*Om2^2)^2-4*Om1^2*Om2^2)));
        lambda1 = om1^2 / sqrt((om1^2-Om2^2)^2+mu*Om2^4);
        lambda2 = om2^2 / sqrt((om2^2-Om2^2)^2+mu*Om2^4);
        D1 = lambda1;
        D2 = lambda2;
        
        G1 = 1/2*D1*lambda1^alpha*om1^(alpha-1)*pi^(-3/2)*gamma(alpha/2+1)/gamma(alpha/2+3/2);
        G2 = 1/2*D2*lambda2^alpha*om2^(alpha-1)*pi^(-3/2)*gamma(alpha/2+1)/gamma(alpha/2+3/2);
        
        
        D = @(A) -epsilon*mu * [G1 * Ialpha(lambda2*om2*A(2)/(lambda1*om1*A(1)), alpha) * A(1)^alpha ;...
            G2 * Ialpha(lambda1*om1*A(1)/(lambda2*om2*A(2)), alpha) * A(2)^alpha];
        
        DA = nan(size(A));
        for i = 1:size(A,2)
            DA(:,i) = D(A(:,i));
        end
    end




%%

    function I = Ialpha(lambda, alpha)
        if alpha<1 && abs(lambda) == inf
            I = 0;
            return
        end
        theta = linspace(0, 2*pi, nIalpha);
        theta = theta(1:end-1);
        I = (1 + lambda^2 + 2*lambda*cos(theta)).^(alpha/2-1/2) .* (1 + lambda*cos(theta));
        I = sum(I)*(theta(2)-theta(1));

%         I = 2*pi * (1+lambda.^2).^((alpha-1)/2) * (1 + (alpha-1)/2*lambda^2/(1+lambda^2));
        
%         I = 2*pi;
    end



%% integration temporelle


nT = 100;

Xi = [xi(1); xi(1)+xi(2)];
Vi = [vi(1); vi(1)+vi(2)];

Dtemp = @(t, Y) [
    Y(3);
    Y(4);
    -w0^2*Y(1) + mu*w1^2*(Y(2)-Y(1)) + mu*epsilon*sign(Y(4)-Y(3))*abs(Y(4)-Y(3))^alpha;
    -w1^2*(Y(2)-Y(1))-epsilon*sign(Y(4)-Y(3))*abs(Y(4)-Y(3))^alpha
    ];


tnext = 0;
    function status = outputWait2(t, ~, ~)
        status = 0;
        if t>tnext
            waitbar(t(1)/T, wait, ['integration (' num2str(round(t(1)/T*100)) '%)']);
            tnext = tnext + T/200;
        end
    end
    

wait = waitbar(0, 'integration (0%)');

options = odeset('RelTol', 1e-10, 'Stats', 'off', 'MaxStep', 1/(w0*nT), 'OutputFcn', @outputWait2);

[t0, Y] = ode45(Dtemp, [0 T], [Xi; Vi], options);


xtemp = Y(:,1);

t2 = linspace(t0(1), t0(end), T*w0/(2*pi)*nT);
xtemp = interp1(t0, xtemp, t2);

t0 = t2;


close(wait);


fig = figure;
ax = axes(fig);
waveletplot = plot(t0, xtemp, 'Parent', ax);
grid(ax, 'on');
xlabel(ax, '$t$');
ylabel(ax, '$x_1$');





%% ondelette


Q = 25;
MaxRidges = 2;
MaxParallelRidges = 2;
fmin = 0.9;
fmax = 1.1;
NbFreq = 100;


WaveletMenu('WaveletPlot', waveletplot, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxRidges', MaxRidges, 'MaxParallelRidges', MaxParallelRidges);


%% extraction de dA
nReg = 100;

RegressionMenu;

[t01, At1] = PlotExtract;
[t02, At2] = PlotExtract;

dAt1 = zeros(size(At1));
for k = 1:length(At1)-1
    dAt1(k) = (At1(k+1) - At1(k))/(t01(k+1) - t01(k));
end
for k = 2:length(At1)
    dAt1(k) = dAt1(k) + (At1(k) - At1(k-1))/(t01(k) - t01(k-1));
end
dAt1(2:end-1) = dAt1(2:end-1)/2;

dAt2 = zeros(size(At2));
for k = 1:length(At2)-1
    dAt2(k) = (At2(k+1) - At2(k))/(t02(k+1) - t02(k));
end
for k = 2:length(At2)
    dAt2(k) = dAt2(k) + (At2(k) - At2(k-1))/(t02(k) - t02(k-1));
end
dAt2(2:end-1) = dAt2(2:end-1)/2;

t0 = linspace(max(t01(1), t02(1)), min(t01(end), t02(end)), nReg);

dAt1 = interp1(t01, dAt1, t0);
dAt2 = interp1(t02, dAt2, t0);

dAt = [dAt1; dAt2];

At1 = interp1(t01, At1, t0);
At2 = interp1(t02, At2, t0);

At = [At1; At2];

%% inversion
Om1 = w0;
Om2 = w1;
om1 = sqrt(1/2 * (Om1^2+(1+mu)*Om2^2 - sqrt((Om1^2+(1+mu)*Om2^2)^2-4*Om1^2*Om2^2)));
om2 = sqrt(1/2 * (Om1^2+(1+mu)*Om2^2 + sqrt((Om1^2+(1+mu)*Om2^2)^2-4*Om1^2*Om2^2)));



Mxz = abs(diag([sqrt((Om2^2-om1^2)^2+mu*Om2^4)/(Om2^2-om1^2), sqrt((Om2^2-om2^2)^2+mu*Om2^4)/(Om2^2-om2^2)]));
At = Mxz*At;
dAt = Mxz*dAt;



%%

om1t = 2*pi * 0.94655;
om2t = 2*pi * 1.0459;

%% nb de points regression
%nPointsReg = 1000;

%%
c1 = om1t^2 + om2t^2;
c2 = (om1t^2 + om2t^2)^2 - (om1t^2 - om2t^2)^2;



%% regression

    function deltax = S(Param)
        Pmu = 0.01;
        Pepsilon = Param(1);
        Palpha = Param(2);
        
        Om2t = sqrt(1/2/(1+Pmu) * (c1 - sqrt(c1^2-(1+Pmu)*c2)));
        Om1t = sqrt(c2)/2/Om2t;
        
        DA = derivAprox(mu, w0, w1, Pepsilon, Palpha, At);
        deltax = DA-dAt;
        deltax = deltax(:);
    end


%--test--
param1 = linspace(0, 0.3, 30);
param2 = linspace(0, 2, 30);
[param13D, param23D] = meshgrid(param1, param2);

Z = nan(length(param1), length(param2));
for ind1 = 1:length(param1)
    for ind2 = 1:length(param2)
        Z(ind1,ind2) = sum (S([param1(ind1), param2(ind2)]).^2);
    end
end

figure;
surf(param13D, param23D, log(Z));





%--test--



Param0 = [0.1, 2]; %epsilon, alpha

lb = [0, 0];
ub = [2, 2];

options = optimoptions(@lsqnonlin, 'FunctionTolerance', 0, 'MaxFunctionEvaluations', 1e4,...
    'MaxIterations', 1e4, 'FiniteDifferenceStepSize', 1e-6);

Param1 = lsqnonlin(@S, Param0, lb, ub, options)
















end

