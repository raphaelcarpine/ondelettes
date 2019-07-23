function polesNonLin()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

precisionNumerique = 18;
printReponseTemp = true;

nIalpha = 1000;

%%

mu = 0.01;
w0 = 2*pi;
w1 = 2*pi;
epsilon = 3;
alpha = 1.5;


mu = 0.01;
w0 = 2*pi;
w1 = 2*pi;
epsilon = 0.05;
alpha = 1.9;



T = 100;
dt = 1e-2;

x0 = [0; 0];
v0 = [1; -1/sqrt(mu)];


%% conditions initiales

WaWb = InitPolesNonLin(mu, w0, w1, epsilon, alpha, x0, v0);
wa = WaWb(1);
wb = WaWb(2);

deformA0 = -(w0^2+(1+mu)*wa^2)/(mu*wa^2);
deformB0 = -(w0^2+(1+mu)*wb^2)/(mu*wb^2);

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

% for iterations = 1:1000 % initialisation des CI (on part des valeurs non amorties puis on itere avec l'amortissement)
%     dphidpsi = dPhidPsi([phi0, psi0], [wa, wb], [deformA0, deformB0]);
%     wa = dphidpsi(1);
%     wb = dphidpsi(2);
%     deformA0 = -(w0^2+(1+mu)*wa^2)/(mu*wa^2);
%     deformB0 = -(w0^2+(1+mu)*wb^2)/(mu*wb^2);
%     [phi0, psi0] = Phi0Psi0(wa, wb, deformA0, deformB0);
% end


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
            * Ialpha(abs(imag(dpsi)/imag(dphi))*exp(real(psi-phi))*abs(deformB/deformA));
        Cb = epsilon*1i/pi*exp(real(psi)*(alpha-1))*imag(dpsi)^alpha*abs(deformB)^(alpha-1)...
            * gamma(alpha/2+1)/(sqrt(pi)*gamma(alpha/2+3/2))...
            * Ialpha(abs(imag(dphi)/imag(dpsi))*exp(real(phi-psi))*abs(deformA/deformB));
        if isnan(Ca)
            Ca = 0;
        end
        if isnan(Cb)
            Cb = 0;
        end
        
        
        dphi2 = -(w0^2+(1+mu)*w1^2+(1+mu)*Ca)/2 - 1/2*sqrt((w0^2+(1+mu)*w1^2+(1+mu)*Ca)^2-4*w0^2*(w1^2+Ca));
        dpsi2 = -(w0^2+(1+mu)*w1^2+(1+mu)*Cb)/2 + 1/2*sqrt((w0^2+(1+mu)*w1^2+(1+mu)*Cb)^2-4*w0^2*(w1^2+Cb));
        
        
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


%% integration

tout = 0:dt:T;
dAnglesout = nan(length(tout), 2);
dAnglesout(1,:) = [wa, wb];
Anglesout = nan(length(tout), 2);
Anglesout(1,:) = [phi0, psi0];
Deforms = nan(length(tout), 2);
Deforms(1,:) = [deformA0, deformB0];

wait = waitbar(0, 'integration 1/2 (0%)');
for k = 2:length(tout)
    if k == 2
        dAnglesout(k-1,:) = dPhidPsi(Anglesout(k-1,:), dAnglesout(k-1,:), Deforms(k-1,:));
    end
    Anglesout(k,:) = Anglesout(k-1,:) + dt*dAnglesout(k-1,:);
    Deforms(k,:) = - (w0^2 + (1+mu)*dAnglesout(k-1,:).^2) ./ (mu*dAnglesout(k-1,:).^2);
    
    dAnglesout(k,:) = dPhidPsi(Anglesout(k,:), dAnglesout(k-1,:), Deforms(k,:));
    Anglesout(k,:) = Anglesout(k-1,:) + dt*(dAnglesout(k-1,:)+dAnglesout(k,:))/2;
    
    if real(Anglesout(k,1)) < real(Anglesout(k,2)) - precisionNumerique
        Anglesout(k,1) = -inf + 1i*imag(Anglesout(k,1));
    end
    if real(Anglesout(k,2)) < real(Anglesout(k,1)) - precisionNumerique
        Anglesout(k,2) = -inf + 1i*imag(Anglesout(k,2));
    end
    
    Deforms(k,:) = - (w0^2 + (1+mu)*dAnglesout(k,:).^2) ./ (mu*dAnglesout(k,:).^2);
    
    if mod(k, round(length(tout)/200)) == 0
        waitbar(k/length(tout), wait, ['integration 1/2 (' num2str(round(k/length(tout)*100)) '%)'])
    end
end
% close(wait);


%% integration temporelle


nT = 2000;

X0 = [x0(1); x0(1)+x0(2)];
V0 = [v0(1); v0(1)+v0(2)];

D = @(t, Y) [
    Y(3);
    Y(4);
    -w0^2*Y(1) + mu*w1^2*(Y(2)-Y(1)) + mu*epsilon*sign(Y(4)-Y(3))*abs(Y(4)-Y(3))^alpha;
    -w1^2*(Y(2)-Y(1))-epsilon*sign(Y(4)-Y(3))*abs(Y(4)-Y(3))^alpha
    ];


tnext = 0;
    function status = outputWait(t, ~, ~)
        status = 0;
        if t>tnext
            waitbar(t(1)/T, wait, ['integration 2/2 (' num2str(round(t(1)/T*100)) '%)']);
            tnext = tnext + T/200;
        end
    end
    
if printReponseTemp
    
    waitbar(0, wait, 'integration 2/2 (0%)');
    
    options = odeset('RelTol', 1e-10, 'Stats', 'off', 'MaxStep', 1/(w0*nT), 'OutputFcn', @outputWait);
    
    [t, Y] = ode45(D, [0 T], [X0; V0], options);
    
    
    X0 = Y(:,1);
    X1 = Y(:,2);
    
    t2 = linspace(t(1), t(end), T*nT);
    X0 = interp1(t, X0, t2);
    X1 = interp1(t, X1, t2);
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
if printReponseTemp
    plot(t, X1-X0, 'Parent', ax);
end
plot(tout, real(sum(exp(Anglesout).*Deforms, 2)), 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylabel(ax, 'x1');

fig = figure;
ax = axes(fig);
hold(ax, 'on');
if printReponseTemp
    waveletplot = plot(t, X0, 'Parent', ax);
end
% plot(t2, X2, 'Parent', ax);
plot(tout, real(sum(exp(Anglesout), 2)), 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylabel(ax, 'x0');
% ylim(ax, [-1, 1]);

fig = figure;
ax = axes(fig);
plot(tout, imag(dAnglesout)/(2*pi), 'Parent', ax);
grid(ax, 'on');
ylabel(ax, 'freqs');
% ylim(ax, [0, 2]);

fig = figure;
ax = axes(fig);
plot(tout, -real(dAnglesout), 'Parent', ax);
grid(ax, 'on');
ylabel(ax, '\lambda');
% ylim(ax, [-1, 1]);

fig = figure;
ax = axes(fig);
plot(tout, exp(real(Anglesout)), tout, sum(exp((1-alpha)*real(Anglesout)), 2), 'Parent', ax);
grid(ax, 'on');
ylabel(ax, 'abs ridges x0');

fig = figure;
ax = axes(fig);
plot(tout, abs(exp(Anglesout).*Deforms), tout, sum(abs(exp((1-alpha)*Anglesout).*Deforms), 2), 'Parent', ax);
grid(ax, 'on');
ylabel(ax, 'abs ridges x1');

% fig = figure;
% ax = axes(fig);
% plot(tout, angle(Deforms), 'Parent', ax);
% grid(ax, 'on');
% ylabel(ax, 'Im(\delta)');
% 
% fig = figure;
% ax = axes(fig);
% plot(tout, abs(Deforms), 'Parent', ax);
% grid(ax, 'on');
% ylabel(ax, '|e^{\delta}|');

diffAngles = abs(imag(dAnglesout(:,1)) - imag(dAnglesout(:,2)));
ddiff = abs([diff(diffAngles); 0]./diffAngles)/dt;
dDeforma = abs([diff(Deforms(:,1)); 0]./Deforms(:,1))/dt;
dDeformb = abs([diff(Deforms(:,2)); 0]./Deforms(:,2))/dt;
ddAnglesa = abs([diff(dAnglesout(:,1)); 0]./dAnglesout(:,1))/dt;
ddAnglesb = abs([diff(dAnglesout(:,2)); 0]./dAnglesout(:,2))/dt;


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(tout, diffAngles, 'Parent', ax);
plot(tout, ddiff, 'Parent', ax);
plot(tout, dDeforma, 'Parent', ax);
plot(tout, dDeformb, 'Parent', ax);
plot(tout, ddAnglesa, 'Parent', ax);
plot(tout, ddAnglesb, 'Parent', ax);
hold(ax, 'off');
grid(ax, 'on');
ylabel(ax, '');
legend('|\omega - \Omega|', "|(\omega - \Omega)'/(\omega - \Omega)|", "|\delta'|", "|\Delta'|",...
    "|\omega'/\omega|", "|\Omega'/\Omega|");


%% ondelette


Q = 25;
MaxParallelRidges = 2;
fmin = 0.9;
fmax = 1.1;
NbFreq = 100;


WaveletMenu('WaveletPlot', waveletplot, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);


end

