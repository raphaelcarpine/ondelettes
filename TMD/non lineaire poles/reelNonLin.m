function reelNonLin()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

printReponseTemp = true;

nIalpha = 1000;

%%

mu = 0.01;
w0 = 2*pi;
w1 = 2*pi/(1+mu);
epsilon = 2*w1*sqrt(1*mu/1/(1+mu));
alpha = 1.;



mu = 0.01;
w0 = 2*pi;
w1 = 2*pi;
epsilon = 0.1;
alpha = 0.5;



T = 100;
dt = 1e-2;

x0 = [0; 0];
v0 = [0; 1];

%% approx sans les poles

nT = 5;

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

delta10 = epsilon * (om1+om2)/(om2-om1)/2
delta2 = 2*(om2-om1)/(om1+om2)

matM12 = diag([1 sqrt(mu)]);
matO = [(Om2^2-om1^2)/sqrt((om1^2-Om2^2)^2+mu*Om2^4), (Om2^2-om2^2)/sqrt((om2^2-Om2^2)^2+mu*Om2^4);...
    sqrt(mu)*Om2^2/sqrt((om1^2-Om2^2)^2+mu*Om2^4), sqrt(mu)*Om2^2/sqrt((om2^2-Om2^2)^2+mu*Om2^4)];

% conditions initiales
X0 = [x0(1); x0(1)+x0(2)];
V0 = [v0(1); v0(1)+v0(2)];

Za0 = 1i*transpose(matO)*matM12*X0 + diag([om1 om2])\transpose(matO)*matM12*V0;
A0 = abs(Za0);
Phi0 = angle(Za0);

% integration
D = @(t, A) -epsilon*mu * [G1 * Ialpha(lambda2*om2*A(2)/(lambda1*om1*A(1))) * A(1)^alpha ;...
    G2 * Ialpha(lambda1*om1*A(1)/(lambda2*om2*A(2))) * A(2)^alpha];

wait = waitbar(0, 'integration 1/2 (0%)');
tnext = 0; % variable d'affichage
    function status = outputWait(t, ~, ~)
        status = 0;
        if t>tnext
            waitbar(t(1)/T, wait, ['integration 1/2 (' num2str(round(t(1)/T*100)) '%)']);
            tnext = tnext + T/200;
        end
    end

options = odeset('RelTol', 1e-10, 'Stats', 'off', 'OutputFcn', @outputWait, 'MaxStep', 1/(w0*nT));

[tr, A] = ode45(D, [0 T], A0, options);

A = transpose(A);
tr = transpose(tr);

Zr = A .* sin([om1*tr+Phi0(1); om2*tr+Phi0(2)]);
Xr = matM12\matO*Zr;
X0r = Xr(1,:);
X1r = Xr(2,:);



%%

    function I = Ialpha(lambda)
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


nT = 2000;

X0 = [x0(1); x0(1)+x0(2)];
V0 = [v0(1); v0(1)+v0(2)];

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
            waitbar(t(1)/T, wait, ['integration 2/2 (' num2str(round(t(1)/T*100)) '%)']);
            tnext = tnext + T/200;
        end
    end
    
if printReponseTemp
    
    waitbar(0, wait, 'integration /2/ (0%)');
    
    options = odeset('RelTol', 1e-10, 'Stats', 'off', 'MaxStep', 1/(w0*nT), 'OutputFcn', @outputWait2);
    
    [t, Y] = ode45(Dtemp, [0 T], [X0; V0], options);
    
    
    X0 = Y(:,1);
    X1 = Y(:,2);
    V0 = Y(:,3);
    V1 = Y(:,4);
    
    t2 = linspace(t(1), t(end), T*nT);
    X0 = interp1(t, X0, t2);
    X1 = interp1(t, X1, t2);
    V0 = interp1(t, V0, t2);
    V1 = interp1(t, V1, t2);
    Xtemp = [X0; X1];
    Vtemp = [V0; V1];
    Atemp = abs (1i*transpose(matO)*matM12*Xtemp + diag([om1 om2])\transpose(matO)*matM12*Vtemp);
    t = t2;
end


close(wait);




%% integration temporelle 2
% nT = 1000;
% 
% 
% delta = 0.01;
% lambda = epsilon/delta/w0*beta(1/2, alpha/2+1)/2;
% 
% D2 = @(t, Y) [
%     Y(3);
%     Y(4);
%     -w0^2*Y(1) + mu*w1^2*(Y(2)-Y(1)) + mu*lambda*(abs(Y(2)-Y(1))<=delta)*(Y(4)-Y(3));
%     -w1^2*(Y(2)-Y(1))-lambda*(abs(Y(2)-Y(1))<=delta)*(Y(4)-Y(3))
%     ];
% 
% 
% waitbar(0, wait, 'integration 2/2 (0%)');
% 
% tnext = 0;
% options = odeset('RelTol', 1e-10, 'Stats', 'off', 'MaxStep', 1/(w0*nT), 'OutputFcn', @outputWait);
% 
% [t2, Y2] = ode45(D2, [0 T], [X0; V0], options);
% 
% 
% X2 = Y2(:,1);
% 
% close(wait);



%% affichage

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(t, X1-X0, '--', 'Parent', ax);
plot(tr, X1r-X0r, 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylabel(ax, 'x1');

fig = figure;
ax = axes(fig);
hold(ax, 'on');
waveletplot = plot(t, X0, '--', 'Parent', ax);
plot(tr, X0r, 'Parent', ax);
% plot(t2, X2, 'Parent', ax);
%plot(tout, real(sum(exp(Anglesout), 2)), 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylabel(ax, 'x0');
% ylim(ax, [-1, 1]);

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(t, Atemp(1,:), '--', 'Parent', ax);
plot(t, Atemp(2,:), '--', 'Parent', ax);
plot(tr, A(1,:), 'Parent', ax);
plot(tr, A(2,:), 'Parent', ax);
grid(ax, 'on');
ylabel(ax, 'abs ridges');

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(t, Atemp(1,:).^(1-alpha) + Atemp(2,:).^(1-alpha), '--', 'Parent', ax);
plot(tr, A(1,:).^(1-alpha) + A(2,:).^(1-alpha), 'Parent', ax);
grid(ax, 'on');
ylabel(ax, 'sum A^(1-alpha)');


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(t, sum(Atemp .* (diag([om1^2, om2^2]) * Atemp), 1), '--', 'Parent', ax);
plot(tr, sum(A .* (diag([om1^2, om2^2]) * A), 1), 'Parent', ax);
grid(ax, 'on');
ylabel(ax, 'énergie');


%% ondelette


Q = 25;
MaxParallelRidges = 2;
fmin = 0.9;
fmax = 1.1;
NbFreq = 100;


WaveletMenu('WaveletPlot', waveletplot, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);


end

