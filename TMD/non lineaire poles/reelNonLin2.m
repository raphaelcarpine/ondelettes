function reelNonLin2()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

printReponseTemp = true;


%%
mu = 0.01;
w0 = 2*pi;
w1 = 2*pi/(1+mu);
epsilon = 1*w1*sqrt(3*mu/8/(1+mu));
alpha = 1.;



% mu = 0.01;
% w0 = 2*pi;
% w1 = 2*pi/(1+mu);
% epsilon = 10.5;
% alpha = 0.5;



T = 50;
dt = 1e-2;

x0 = [0; 0];
v0 = [1; 0];

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
lw1 = lambda1*om1;
lw2 = lambda2*om2;

delta10 = epsilon
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
D = @(t, Aphi) - epsilon * mu...
    * ((lw1*Aphi(1))^2+(lw2*Aphi(2))^2+2*lw1*Aphi(1)*lw2*Aphi(2)*cos((om2-om1)*t+Aphi(4)-Aphi(3)))^((alpha-1)/2)...
    * gamma(alpha/2+1) / (sqrt(pi) * gamma(alpha/2+3/2)) * [
    D1/om1 * (lw1*Aphi(1) + lw2*Aphi(2) * cos((om2-om1)*t+Aphi(4)-Aphi(3)));
    D2/om2 * (lw2*Aphi(2) + lw1*Aphi(1) * cos((om2-om1)*t+Aphi(4)-Aphi(3)));
     D1/(om1*Aphi(1)) * lw2*Aphi(2) * sin((om2-om1)*t+Aphi(4)-Aphi(3));
    -D2/(om2*Aphi(2)) * lw1*Aphi(1) * sin((om2-om1)*t+Aphi(4)-Aphi(3))
    ];


wait = waitbar(0, 'integration 1/2 (0%)');
tnext = 0; % variable d'affichage
    function status = outputWait(t, ~, ~)
        status = 0;
        if t>tnext
            waitbar(t(1)/T, wait, ['integration 1/2 (' num2str(round(t(1)/T*100)) '%)']);
            tnext = tnext + T/200;
        end
    end

options = odeset('RelTol', 1e-10, 'OutputFcn', @outputWait, 'MaxStep', 1/(w0*nT));

[tr, Aphi] = ode45(D, [0 T], [A0; Phi0], options);

Aphi = transpose(Aphi);
tr = transpose(tr);

A = Aphi(1:2,:);
Phi = Aphi(3:4,:);
dPhi = nan(size(Aphi));
for i = 1:size(Aphi,2)
    dPhi(:,i) = D(tr(i), Aphi(:,i));
end
dPhi = dPhi(3:4,:);

Zr = A .* sin([om1*tr+Phi(1,:); om2*tr+Phi(2,:)]);
Xr = matM12\matO*Zr;
X0r = Xr(1,:);
X1r = Xr(2,:);





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
    
    options = odeset('RelTol', 1e-10, 'MaxStep', 1/(w0*nT), 'OutputFcn', @outputWait2);
    
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
    t = t2;
    
    
    Xtemp = [X0; X1];
    Vtemp = [V0; V1];
    Ztemp = 1i*transpose(matO)*matM12*Xtemp + diag([om1 om2])\transpose(matO)*matM12*Vtemp;
    Atemp = abs (Ztemp);
    Phitemp = angle (Ztemp .* exp(-1i * diag([om1 om2]) * [t; t]));
    for k = 1:size(Phitemp, 2)-1
        for l = 1:2
            if round((Phitemp(l,k+1)-Phitemp(l,k))/(2*pi)) ~= 0
                Phitemp(l,k+1:end) = Phitemp(l,k+1:end) - 2*pi * round((Phitemp(l,k+1)-Phitemp(l,k))/(2*pi));
            end
        end
    end
    dPhitemp = diff(Phitemp, 1, 2);
    dPhitemp = [dPhitemp(:,1), 1/2*(dPhitemp(:,1:end-1)+dPhitemp(:,2:end)), dPhitemp(:,end)];
    dPhitemp = dPhitemp/(t(2)-t(1));
end


close(wait);

%% ridge unique

matZX = matM12\matO;

cX0r = A .* transpose(matZX(1,:)); % composantes de X0r
cX1r = A .* transpose(matZX(2,:));
ridgeX0r = sqrt(cX0r(1,:).^2 + cX0r(2,:).^2 + 2*cX0r(1,:).*cX0r(2,:).*cos((om2-om1)*tr+Phi(2,:)-Phi(1,:)));
ridgeX1r = sqrt(cX1r(1,:).^2 + cX1r(2,:).^2 + 2*cX1r(1,:).*cX1r(2,:).*cos((om2-om1)*tr+Phi(2,:)-Phi(1,:)));

if printReponseTemp
    cX0t = Atemp .* transpose(matZX(1,:)); % composantes de X0
    cX1t = Atemp .* transpose(matZX(2,:));
    ridgeX0t = sqrt(cX0t(1,:).^2 + cX0t(2,:).^2 + 2*cX0t(1,:).*cX0t(2,:).*cos((om2-om1)*t+Phitemp(2,:)-Phitemp(1,:)));
    ridgeX1t = sqrt(cX1t(1,:).^2 + cX1t(2,:).^2 + 2*cX1t(1,:).*cX1t(2,:).*cos((om2-om1)*t+Phitemp(2,:)-Phitemp(1,:)));
end





%% affichage


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tr, X1r-X0r, 'Parent', ax);
if printReponseTemp
    plot(t, X1-X0, 'Parent', ax);
end
%plot(tout, real(sum(exp(Anglesout).*Deforms, 2)), 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylabel(ax, 'x1');

fig = figure;
ax = axes(fig);
hold(ax, 'on');
waveletplot = plot(tr, X0r, 'b', 'Parent', ax);
if printReponseTemp
    waveletplot = plot(t, X0, 'r', 'Parent', ax);
    uistack(waveletplot, 'bottom');
end
% plot(t2, X2, 'Parent', ax);
%plot(tout, real(sum(exp(Anglesout), 2)), 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
xlabel(ax, '$t$');
ylabel(ax, '$x_1$');
% ylim(ax, [-1, 1]);

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tr, A(1,:), 'b', 'Parent', ax);
plot(tr, A(2,:), 'b', 'Parent', ax);
if printReponseTemp
    uistack(plot(t, Atemp(1,:), 'r', 'Parent', ax), 'bottom');
    uistack(plot(t, Atemp(2,:), 'r', 'Parent', ax), 'bottom');
end
grid(ax, 'on');
ylabel(ax, 'abs ridges');


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tr, ridgeX0r, 'Parent', ax);
plot(tr, X0r, 'Parent', ax);
if printReponseTemp
    plot(t, ridgeX0t, 'Parent', ax);
    plot(t, X0, 'Parent', ax);
end
grid(ax, 'on');
ylabel(ax, 'abs ridge unique');


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tr, lw1*A(1,:), 'Parent', ax);
plot(tr, lw2*A(2,:), 'Parent', ax);
plot(tr, lw2*A(2,:) - lw1*A(1,:), 'Parent', ax);
grid(ax, 'on');
ylabel(ax, 'a2-a1');


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tr, Phi(1,:), 'Parent', ax);
plot(tr, Phi(2,:), 'Parent', ax);
plot(tr, Phi(2,:)+Phi(1,:), 'Parent', ax);
if printReponseTemp
    plot(t, Phitemp(1,:), 'Parent', ax);
    plot(t, Phitemp(2,:), 'Parent', ax);
    plot(t, Phitemp(2,:)+Phitemp(1,:), 'Parent', ax);
end
grid(ax, 'on');
ylabel(ax, '\phi_1, \phi_2, \phi_1+\phi_2');


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tr, (om2-om1)*tr + Phi(2,:) - Phi(1,:), 'Parent', ax);
plot(tr, cos((om2-om1)*tr + Phi(2,:) - Phi(1,:)), 'Parent', ax);
plot(tr, sin((om2-om1)*tr + Phi(2,:) - Phi(1,:)), 'Parent', ax);
if printReponseTemp
    plot(t, (om2-om1)*t + Phitemp(2,:) - Phitemp(1,:), 'Parent', ax);
    plot(t, cos((om2-om1)*t + Phitemp(2,:) - Phitemp(1,:)), 'Parent', ax);
    plot(t, sin((om2-om1)*t + Phitemp(2,:) - Phitemp(1,:)), 'Parent', ax);
end
grid(ax, 'on');
ylabel(ax, '\phi, \cos\phi, \sin\phi');


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tr, dPhi(1,:) + om1, 'Parent', ax);
plot(tr, dPhi(2,:) + om2, 'Parent', ax);
if printReponseTemp
    plot(t, dPhitemp(1,:) + om1, 'Parent', ax);
    plot(t, dPhitemp(2,:) + om2, 'Parent', ax);
end
grid(ax, 'on');
ylabel(ax, 'freq instantanée');


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tr, sum(A .* (diag([om1^2, om2^2]) * A), 1), 'Parent', ax);
if printReponseTemp
    plot(t, sum(Atemp .* (diag([om1^2, om2^2]) * Atemp), 1), 'Parent', ax);
end
grid(ax, 'on');
ylabel(ax, 'énergie');

% fig = figure;
% ax = axes(fig);
% hold(ax, 'on');
% plot(tr, A(1,:).^(1-alpha) + A(2,:).^(1-alpha), 'Parent', ax);
% if printReponseTemp
%     plot(t, Atemp(1,:).^(1-alpha) + Atemp(2,:).^(1-alpha), 'Parent', ax);
% end
% grid(ax, 'on');
% ylabel(ax, 'sum A^(1-alpha)');

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tr, (A(1,:)-A(2,:)).^2, 'Parent', ax);
plot(tr, A(1,:).^2+A(2,:).^2+2*A(1,:).*A(2,:).*cos((om2-om1)*tr + Phi(2,:)-Phi(1,:)) - (A(1,:)-A(2,:)).^2, 'Parent', ax);
plot(tr, A(1,:)-A(2,:), 'Parent', ax);
plot(tr, A(2,:)+A(2,:).*cos((om2-om1)*tr + Phi(2,:)-Phi(1,:)), 'Parent', ax);
grid(ax, 'on');
legend();
ylabel(ax, '(a_2-a_1)^2, a_1 a_2 \delta^2');


%% ondelette


Q = 25;
MaxParallelRidges = 2;
fmin = 0.9;
fmax = 1.1;
NbFreq = 100;


WaveletMenu('WaveletPlot', waveletplot, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);


end

