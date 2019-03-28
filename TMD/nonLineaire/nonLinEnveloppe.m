function nonLinEnveloppe()
clear all;
% close all;



mu = 0.01;
omega0 = 2*pi;
zeta0 = 0.0;
omega1 = omega0/(1+mu);
omega1 = omega0;
a = @(x1, v1) 0.03*2*omega1*sign(v1)*abs(v1)^0.5;

T = 100;
nbPointsPeriode = 100;
nT = T*nbPointsPeriode*omega0/2/pi;

x0 = [0; 0];
v0 = [1; -1];

%%
w0 = omega0;
w02 = omega0^2;
w1 = omega1;
w12 = omega1^2;
z0 = zeta0;

%%
t = linspace(0, T, 10000);

d2x0 = @(x0, x1, dx0, dx1) -2*z0*w0*dx0 - w02*x0 + mu*a(x1,dx1) + mu*w12*x1;
d2x1 = @(x0, x1, dx0, dx1) -(1+mu)*a(x1,dx1) - (1+mu)*w12*x1 + 2*z0*w0*dx0 + w02*x0;

DX = @(t, XdX) [XdX(3); XdX(4); d2x0(XdX(1),XdX(2),XdX(3),XdX(4)) ; d2x1(XdX(1),XdX(2),XdX(3),XdX(4))];

tic;
[tout, Xout] = differentialEq([x0; v0], DX, T, false);
toc

t = linspace(0, T, nT);
X = interp1(tout, Xout, t);


f = figure;
ax = axes(f);
ax.XGrid = 'on';
ax.YGrid = 'on';
hold(ax, 'on');

% reponseTemp = [plot(t, X(:,1), 'Parent', ax), plot(t, sqrt(mu)*X(:,2), 'Parent', ax)];
reponseTemp = plot(t, X(:,1), 'Parent', ax);

tic;
hilb = hilbert([X(:,1); zeros(length(X(:,1)), 1)]);
hilb = hilb(1:length(X(:,1)));
toc
plot(t, -abs(hilb), 'Parent', ax);


%%
A = @(phi1,dphi1) a(real(exp(1i*phi1)), real(1i*dphi1*exp(1i*phi1)))...
    + 1i*a(imag(exp(1i*phi1)), imag(1i*dphi1*exp(1i*phi1)));

d2phi0 = @(p0,p1,dp0,dp1) 1i*( -dp0^2 + 2*z0*w0*1i*dp0 + w02 - mu*A(p1,dp1)*exp(-1i*p0) - mu*w12*exp(1i*(p1-p0)));
d2phi1 = @(p0,p1,dp0,dp1) 1i*( -dp1^2 + (1+mu)*A(p1,dp1)*exp(-1i*p1) + (1+mu)*w12...
    -(2*z0*w0*1i*dp0 + w02)*exp(1i*(p0-p1)));

DPhi = @(t, PdP) [PdP(3); PdP(4); d2phi0(PdP(1),PdP(2),PdP(3),PdP(4)) ; d2phi1(PdP(1),PdP(2),PdP(3),PdP(4))];

Angles0 = [1i*log(omega0)-pi/2 ; 1i*log(omega0)+pi/2; omega1; omega1];


% dangles0 = [v0(1)
% Angles0 = [angles0; dangles0];


tic;
[t, Phi] = differentialEq(Angles0, DPhi, T, false, 'MaxStep', 1e-2);
toc


phi = Phi(:,1);
Xphi = abs(exp(1i*phi));

plot(t, -Xphi, 'Parent', ax);


hold(ax, 'off');


% Rphi = real(phi);
% f2 = figure;
% ax2 = axes(f2);
% plot(t, Rphi, 'Parent', ax2);

%%
Q = 1;
fmin = 0.9;
fmax = 1.1;
NbFreq = 200;




WaveletMenu('fmin',fmin,'fmax',fmax,'NbFreq',NbFreq, 'WaveletPlot', reponseTemp, 'Q', Q);

end

