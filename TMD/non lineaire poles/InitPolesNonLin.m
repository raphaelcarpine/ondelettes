function WaWb = InitPolesNonLin(mu, w0, w1, epsilon, alpha, x0, v0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nIalpha = 1000;
Nmax = 1000;
omegaConvergenceTol = 1e-4;


%% conditions initiales

wa2 = -(w0^2+(1+mu)*w1^2)/2 - 1/2*sqrt(w0^4+(1+mu)^2*w1^4+2*(1+mu)*w0^2*w1^2-4*w0^2*w1^2);
wb2 = -(w0^2+(1+mu)*w1^2)/2 + 1/2*sqrt(w0^4+(1+mu)^2*w1^4+2*(1+mu)*w0^2*w1^2-4*w0^2*w1^2);
wa = sqrt(wa2);
wb = sqrt(wb2);
% deformA0 = -wa2/(wa2+w1^2);
% deformB0 = -wb2/(wb2+w1^2);
deformA0 = -(w0^2+(1+mu)*wa2)/(mu*wa2);
deformB0 = -(w0^2+(1+mu)*wb2)/(mu*wb2);

% phi0 = log((deformB0*v0(1)+v0(2))/(wa*(deformA0-deformB0)));
% psi0 = log((deformA0*v0(1)+v0(2))/(wb*(deformA0-deformB0)));

    function [phi0, psi0] = Phi0Psi0(wa, wb, Da, Db)
        M = [
            1, 0, 1, 0;
            real(Da), -imag(Da), real(Db), -imag(Db);
            real(wa), -imag(wa), real(wb), -imag(wb);
            real(Da*wa), -imag(Da*wa), real(Db*wb), -imag(Db*wb)
            ];
        expphi0psi0 = M\[x0(1); x0(2); v0(1); v0(2)];
        phi0 = log(expphi0psi0(1) + 1i*expphi0psi0(2));
        psi0 = log(expphi0psi0(3) + 1i*expphi0psi0(4));
    end

[phi0, psi0] = Phi0Psi0(wa, wb, deformA0, deformB0);

Wa = nan(1, Nmax);
Wb = nan(1, Nmax);

QteA = nan(1, Nmax);
QteB = nan(1, Nmax);

for iterations = 1:Nmax % initialisation des CI (on part des valeurs non amorties puis on itere avec l'amortissement)
    Wa(iterations) = wa;
    Wb(iterations) = wb;
    wai = wa;
    wbi = wb;
    
    dphidpsi = dPhidPsi([phi0, psi0], [wa, wb], [deformA0, deformB0]);
    wa = dphidpsi(1);
    wb = dphidpsi(2);
    deformA0 = -(w0^2+(1+mu)*wa^2)/(mu*wa^2);
    deformB0 = -(w0^2+(1+mu)*wb^2)/(mu*wb^2);
    
    [phi0, psi0] = Phi0Psi0(wa, wb, deformA0, deformB0);
    
    if abs(wa-wai)<omegaConvergenceTol && abs(wb-wbi)<omegaConvergenceTol
        WaWb = [wa, wb];
        return
    end
    
    QteA(iterations) = (w0^2+(1+mu)*w1^2-mu*Ca)^2-4*w0^2*w1^2;
    QteB(iterations) = (w0^2+(1+mu)*w1^2-mu*Cb)^2-4*w0^2*w1^2;
    
end

WaWb = [nan, nan];

return

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(real(Wa), 'Parent', ax);
plot(real(Wb), 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylabel(ax, '\Re');

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(imag(Wa), 'Parent', ax);
plot(imag(Wb), 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylabel(ax, '\Im');

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(angle(QteA), 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylim(ax, [-pi, pi]);
ylabel(ax, 'A');

fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(angle(QteB), 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylim(ax, [-pi, pi]);
ylabel(ax, 'B');


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(angle(Wa.^2), 'Parent', ax);
plot(angle(Wb.^2), 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylim(ax, [-pi, pi]);
ylabel(ax, 'B');


%%

    function dphidpsi = dPhidPsi(phiPsi, dphidpsi, deforms)
        phi = phiPsi(1);
        psi = phiPsi(2);
        dphi = dphidpsi(1);
        dpsi = dphidpsi(2);
        deformA = deforms(1);
        deformB = deforms(2);
        
        Ca = epsilon*1i/pi*exp(real(phi)*(alpha-1))*imag(dphi)^alpha*abs(deformA)^(alpha-1)...
            * gamma(alpha/2+1)/(sqrt(pi)*gamma(alpha/2+3/2))...
            * Ialpha(abs(imag(dpsi)/imag(dphi)*exp(psi-phi)*deformB/deformA));
        Cb = epsilon*1i/pi*exp(real(psi)*(alpha-1))*imag(dpsi)^alpha*abs(deformB)^(alpha-1)...
            * gamma(alpha/2+1)/(sqrt(pi)*gamma(alpha/2+3/2))...
            * Ialpha(abs(imag(dphi)/imag(dpsi)*exp(phi-psi)*deformA/deformB));
        
        
        dphi2 = -(w0^2+(1+mu)*w1^2+(1+mu)*Ca)/2 - 1/2*sqrt((w0^2+(1+mu)*w1^2+(1+mu)*Ca)^2-4*w0^2*(w1^2+Ca));
        dpsi2 = -(w0^2+(1+mu)*w1^2+(1+mu)*Cb)/2 + 1/2*sqrt((w0^2+(1+mu)*w1^2+(1+mu)*Cb)^2-4*w0^2*(w1^2+Cb));
        
        
        dphi = - sqrt(dphi2);
        dpsi = - sqrt(dpsi2);
%         dphi = (2*(real(dphi)<=0)-1) * dphi;
%         dpsi = (2*(real(dpsi)<=0)-1) * dpsi;
        dphi = (2*(imag(dphi)>=0)-1) * dphi;
        dpsi = (2*(imag(dpsi)>=0)-1) * dpsi;
        
        
        dphidpsi = [dphi, dpsi];
    end

    function dphidpsi = dPhidPsi_bis(phiPsi, dphidpsi, deforms)
        phi = phiPsi(1);
        psi = phiPsi(2);
        dphi = dphidpsi(1);
        dpsi = dphidpsi(2);
        deformA = deforms(1);
        deformB = deforms(2);
        
        Ca = epsilon*1i/pi*exp(real(phi)*(alpha-1))*imag(dphi)^alpha*deformA*abs(deformA)^(alpha-1)...
            * gamma(alpha/2+1)/(sqrt(pi)*gamma(alpha/2+3/2))...
            * Ialpha(abs(imag(dpsi)/imag(dphi)*exp(psi-phi)*deformB/deformA));
        Cb = epsilon*1i/pi*exp(real(psi)*(alpha-1))*imag(dpsi)^alpha*deformB*abs(deformB)^(alpha-1)...
            * gamma(alpha/2+1)/(sqrt(pi)*gamma(alpha/2+3/2))...
            * Ialpha(abs(imag(dphi)/imag(dpsi)*exp(phi-psi)*deformA/deformB));
        
        
        dphi2 = -(w0^2+(1+mu)*w1^2-mu*Ca)/2 - 1/2*sqrt((w0^2+(1+mu)*w1^2-mu*Ca)^2-4*w0^2*w1^2);
        dpsi2 = -(w0^2+(1+mu)*w1^2-mu*Cb)/2 + 1/2*sqrt((w0^2+(1+mu)*w1^2-mu*Cb)^2-4*w0^2*w1^2);
        
        
        
        dphi = - sqrt(dphi2);
        dpsi = - sqrt(dpsi2);
%         dphi = (2*(real(dphi)<=0)-1) * dphi;
%         dpsi = (2*(real(dpsi)<=0)-1) * dpsi;
        dphi = (2*(imag(dphi)>=0)-1) * dphi;
        dpsi = (2*(imag(dpsi)>=0)-1) * dpsi;
        
        
        dphidpsi = [dphi, dpsi];
    end

    function dphidpsi = dPhidPsi_test(phiPsi, dphidpsi, deforms)
        phi = phiPsi(1);
        psi = phiPsi(2);
        dphi = dphidpsi(1);
        dpsi = dphidpsi(2);
        deformA = deforms(1);
        deformB = deforms(2);
        
        Ca = epsilon*1i/pi*exp(real(phi)*(alpha-1))*imag(dphi)^alpha*deformA*abs(deformA)^(alpha-1)...
            * gamma(alpha/2+1)/(sqrt(pi)*gamma(alpha/2+3/2))...
            * Ialpha(abs(imag(dpsi)/imag(dphi)*exp(psi-phi)*deformB/deformA));
        Cb = epsilon*1i/pi*exp(real(psi)*(alpha-1))*imag(dpsi)^alpha*deformB*abs(deformB)^(alpha-1)...
            * gamma(alpha/2+1)/(sqrt(pi)*gamma(alpha/2+3/2))...
            * Ialpha(abs(imag(dphi)/imag(dpsi)*exp(phi-psi)*deformA/deformB));
        
        dphi2m = -(w0^2+(1+mu)*w1^2-mu*Ca)/2 - 1/2*sqrt((w0^2+(1+mu)*w1^2-mu*Ca)^2-4*w0^2*w1^2);
        dphi2p = -(w0^2+(1+mu)*w1^2-mu*Ca)/2 + 1/2*sqrt((w0^2+(1+mu)*w1^2-mu*Ca)^2-4*w0^2*w1^2);
        dpsi2m = -(w0^2+(1+mu)*w1^2-mu*Cb)/2 - 1/2*sqrt((w0^2+(1+mu)*w1^2-mu*Cb)^2-4*w0^2*w1^2);
        dpsi2p = -(w0^2+(1+mu)*w1^2-mu*Cb)/2 + 1/2*sqrt((w0^2+(1+mu)*w1^2-mu*Cb)^2-4*w0^2*w1^2);
        
        
        dphim = - sqrt(dphi2m);
        dphip = - sqrt(dphi2p);
        dpsim = - sqrt(dpsi2m);
        dpsip = - sqrt(dpsi2p);
%         dphi = (2*(real(dphi)<=0)-1) * dphi;
%         dpsi = (2*(real(dpsi)<=0)-1) * dpsi;
        dphim = (2*(imag(dphim)>=0)-1) * dphim;
        dphip = (2*(imag(dphip)>=0)-1) * dphip;
        dpsim = (2*(imag(dpsim)>=0)-1) * dpsim;
        dpsip = (2*(imag(dpsip)>=0)-1) * dpsip;
        
        if abs(dphim-dphi) < abs(dphip-dphi)
            dphi = dphim;
        else
            dphi = dphip;
        end
        if abs(dpsim-dpsi) < abs(dpsip-dpsi)
            dpsi = dpsim;
        else
            dpsi = dpsip;
        end
        
        
        dphidpsi = [dphi, dpsi];
    end


    function I = Ialpha(lambda)        
        theta = linspace(0, 2*pi, nIalpha);
        theta = theta(1:end-1);
        I = (1 + lambda^2 + 2*lambda*cos(theta)).^(alpha/2-1/2) .* (1 + lambda*cos(theta));
        I = sum(I)*(theta(2)-theta(1));

%         I = 2*pi * (1+lambda^2).^(-1/2) * (1 + -1/2*lambda^2/(1+lambda^2));

%         I = 2*pi;
    end



end

