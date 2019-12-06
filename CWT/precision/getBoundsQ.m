function [Qmin, Qmax, Qz] = getBoundsQ(f, Df, Dt, T)
%GETBOUNDSQ Le, Argoul 2004
%   f = (f1+f2)/2
%   Df = f2 - f1
%   Dt = 1/lambda
%   T = total time

cf = 5;
ct = 5;

Qmin = cf * f/(2*Df);

Qmax = 2*pi*f*T / (4*ct*1/2); %mu_psi \simeq 1/2

lambda = 1/Dt;
zeta = lambda / (2*pi*f);
Qz = 1/sqrt(2) / zeta;
end

