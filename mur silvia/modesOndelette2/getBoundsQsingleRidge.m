function [Qmin, Qmax, Qa] = getBoundsQsingleRidge(f1, f2, lambda1, lambda2, f3, T, ct, cf)
%GETBOUNDSQSINGLERIDGE
%   f1 < f2
%   f3 : closest frequency
%   T = total time

if f1 > f2
    f1 = f1 + f2;
    f2 = f1 - f2;
    f1 = f1 - f2;
end

Df = min(abs(f3-f1), abs(f3-f2));
Qmin = cf * f2/(2*Df);

mu_psi = 1/2; % Q >> 1  =>  mu_psi \simeq 1/2
Qmax = 2*pi*f1*T / (4*ct*mu_psi);

Qa = 1/(2*mu_psi*ct) * 2*pi*f1/(Df + lambda1 + lambda2);
end

