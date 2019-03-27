function testnonLin()
clear all;
% close all;



mu = 0.01;
omega0 = 2*pi;
zeta0 = 0.0;
omega1 = omega0/(1+mu);
zeta1 = 0.;

T = 50;

x0 = [0; 0];
v0 = [1; 0];

%%
t = linspace(0, T, 10000);
[x, v, a, xtmd, vtmd, atmd] = reponseTemporelleSystemeLineaire(x0, v0, t, mu, omega0, omega1, zeta0, zeta1);
f = figure;
ax = axes(f);
hold(ax, 'on');
reponseTemp = plot(t, x, 'Parent', ax);

%%

%     function z2 = positive(z)
%         z2 = z.*(2*(real(z)>0)-1);
%         z2 = z;
%     end

fun = @(t, angles) sqrt([omega0^2 - mu*omega1^2*exp(1i*(angles(2)-angles(1)));...
    omega1^2 + mu*omega1^2 - omega0^2*exp(1i*(angles(1)-angles(2)))]);
angles0 = [1i*log(omega0) ; 1i*log(omega1)-pi];


[t, xout] = differentialEq(angles0, fun, T, false, 'MaxStep', 1e-2);

xout = abs(exp(1i*xout));
plot(t, xout, 'Parent', ax);

%%
fun = @(t, angles) [angles(3); angles(4);...
    1i*(omega0^2-angles(3)^2+1i*2*zeta0*omega0*angles(3)-mu*omega1^2*exp(1i*(angles(2)-angles(1))));...
    1i*((1+mu)*omega1^2-angles(4)^2+1i*2*zeta1*omega1*angles(4)-omega0^2*exp(1i*(angles(1)-angles(2))))];

angles0 = [1i*log(omega0) ; 1i*log(omega1)+pi; omega0; omega1*0];


[t, xout] = differentialEq(angles0, fun, T, false, 'MaxStep', 1e-3);

xout = xout(:,1);
xout = abs(exp(1i*xout));
plot(t, xout, 'Parent', ax);


hold(ax, 'off');

%%
Q = 1;
fmin = 0.9;
fmax = 1.1;
NbFreq = 200;




WaveletMenu('fmin',fmin,'fmax',fmax,'NbFreq',NbFreq, 'WaveletPlot', reponseTemp, 'Q', Q);

end

