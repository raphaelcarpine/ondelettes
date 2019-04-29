function polesNonLin()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu = 0.01;
omega0 = '2*pi';
omega1 = '2*pi';
zeta0 = 0;
epsilon = 0.4;
alpha = 0;

T = 100;
nT = 200;

x0 = [0; 0];
v0 = [1; 0];


%%
w0 = eval(omega0);
w1 = eval(omega1);

wa2 = -(w0^2+(1+mu)*w1^2)/2 - 1/2*sqrt(w0^4+(1+mu)^2*w1^4+2*(1+mu)*w0^2*w1^2-4*w0^2*w1^2);
wb2 = -(w0^2+(1+mu)*w1^2)/2 + 1/2*sqrt(w0^4+(1+mu)^2*w1^4+2*(1+mu)*w0^2*w1^2-4*w0^2*w1^2);
wa = sqrt(wa2);
wb = sqrt(wb2);
deformA0 = -wa2/(wa2+w1^2);
deformB0 = -wb2/(wb2+w1^2);


phi0 = log((deformB0*v0(1)+v0(2))/(wa*(deformA0-deformB0)));
psi0 = log((deformA0*v0(1)+v0(2))/(wb*(deformB0-deformA0)));


% deformA = deformA0;
% deformB = deformB0;

%     function dphidpsi = dPhidPsi(t, PhiPsi)
%         phi = PhiPsi(1);
%         psi = PhiPsi(2);
%         
%         dphi2 = -1/(1+deformA) * (w1^2*deformA + epsilon*2*1i/pi^2*exp(-real(phi))*deformA/abs(deformA)*I0(abs(exp(psi-phi)*deformB/deformA)) );
%         dpsi2 = -1/(1+deformB) * (w1^2*deformB + epsilon*2*1i/pi^2*exp(-real(psi))*deformB/abs(deformB)*I0(abs(exp(phi-psi)*deformA/deformB)) );
%         
%         dphi = - sqrt(dphi2);
%         dpsi = - sqrt(dpsi2);
%         dphi = (2*(real(dphi)<=0)-1) * dphi;
%         dpsi = (2*(real(dpsi)<=0)-1) * dpsi;
% %         dphi = (2*(imag(dphi)>=0)-1) * dphi;
% %         dpsi = (2*(imag(dpsi)>=0)-1) * dpsi;
%         
% %         t
% %         dphi
% %         dpsi
%         deformA = - (w0^2+(1+mu)*dphi2)/(mu*dphi2);
%         deformB = - (w0^2+(1+mu)*dpsi2)/(mu*dpsi2);
%         
%         dphidpsi = [dphi; dpsi];
%     end

    function dphidpsi = dPhidPsi(PhiPsi, Deforms)
        phi = PhiPsi(1);
        psi = PhiPsi(2);
        deformA = Deforms(1);
        deformB = Deforms(2);
        
        dphi2 = -1/(1+deformA) * (w1^2*deformA + epsilon*2*1i/pi^2*exp(-real(phi))*deformA/abs(deformA)*I0(abs(exp(psi-phi)*deformB/deformA)) );
        dpsi2 = -1/(1+deformB) * (w1^2*deformB + epsilon*2*1i/pi^2*exp(-real(psi))*deformB/abs(deformB)*I0(abs(exp(phi-psi)*deformA/deformB)) );
        
        dphi = - sqrt(dphi2);
        dpsi = - sqrt(dpsi2);
%         dphi = (2*(real(dphi)<=0)-1) * dphi;
%         dpsi = (2*(real(dpsi)<=0)-1) * dpsi;
        dphi = (2*(imag(dphi)>=0)-1) * dphi;
        dpsi = (2*(imag(dpsi)>=0)-1) * dpsi;
        
%         dphi = dphi - max(real(dphi),0);
%         dpsi = dpsi - max(real(dpsi),0);
        
        dphidpsi = [dphi, dpsi];
    end


    function I = I0(lambda)
        theta = linspace(0, 2*pi, 1000);
        theta = theta(1:end-1);
        I = (1 + lambda^2 + 2*lambda*cos(theta)).^(-1/2) .* (1 + lambda*cos(theta));
        I = sum(I)*(theta(2)-theta(1));
%         I = 2*pi * (1+lambda^2).^(-1/2) * (1 + -1/2*lambda^2/(1+lambda^2));
%         I = 4;
%         I = 0;
    end


% [tout, Anglesout] = ode45(@dPhidPsi, [0 T], [phi0; psi0],...
%     odeset('RelTol', 1e-10, 'Stats', 'off', 'MaxStep', 1/nT));


dt = 1/nT;

tout = 0:dt:T;
dAnglesout = nan(length(tout), 2);
Anglesout = nan(length(tout), 2);
Anglesout(1,:) = [phi0, psi0];
Deforms = nan(length(tout), 2);
Deforms(1,:) = [deformA0, deformB0];

for k = 2:length(tout)
    if any(isnan(dAnglesout(k-1,:)))
        dAnglesout(k-1,:) = dPhidPsi(Anglesout(k-1,:), Deforms(k-1,:));
    end
    Anglesout(k,:) = Anglesout(k-1,:) + dt*dAnglesout(k-1,:);
    Deforms(k,:) = - (w0^2 + (1+mu)*dAnglesout(k-1,:).^2) ./ (mu*dAnglesout(k-1,:).^2);
    
    dAnglesout(k,:) = dPhidPsi(Anglesout(k,:), Deforms(k,:));
    Anglesout(k,:) = Anglesout(k-1,:) + dt*(dAnglesout(k-1,:)+dAnglesout(k,:))/2;
    Deforms(k,:) = - (w0^2 + (1+mu)*dAnglesout(k,:).^2) ./ (mu*dAnglesout(k,:).^2);
end

fig = figure;
ax = axes(fig);
plot(tout, imag(dAnglesout)/(2*pi), 'Parent', ax);
ylim(ax, [0, 2]);

fig = figure;
ax = axes(fig);
plot(tout, real(dAnglesout), 'Parent', ax);
ylim(ax, [-1, 1]);

fig = figure;
ax = axes(fig);
plot(tout, exp(real(Anglesout)), 'Parent', ax);

fig = figure;
ax = axes(fig);
plot(tout, angle(Deforms), 'Parent', ax);

fig = figure;
ax = axes(fig);
plot(tout, abs(Deforms), 'Parent', ax);



















%% integration classique
% d2x0 = '-omega0^2*x0 - 2*omega0*zeta0*dx0 + mu*omega1^2*(x1-x0) + mu*epsilon*sign(dx1-dx0)*abs(dx1-dx0)^alpha';
% d2x1 = '-omega1^2*(x1-x0)-epsilon*sign(dx1-dx0)*abs(dx1-dx0)^alpha';
% 
% waveletplots = systemeQuelconque({'x0', 'x1'}, {d2x0, d2x1}, {'mu', 'omega0', 'omega1', 'zeta0', 'epsilon', 'alpha'},...
%     {mu, omega0, omega1, zeta0, epsilon, alpha},...
%     x0, v0, false, 'T', T, 'nT', nT);
% 
% 
% %ondelette
% Q = 1;
% MaxParallelRidges = 1;
% fmin = 0.9;
% fmax = 1.1;
% NbFreq = 100;
% 
% 
% WaveletMenu('WaveletPlot', waveletplots, 'fmin', fmin, 'fmax', fmax,...
%     'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);
% 
% 
% %regression
% RegressionMenu;



end

