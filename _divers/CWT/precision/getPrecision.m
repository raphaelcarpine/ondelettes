function [Df, Dt] = getPrecision(f, Q)
%GETPRECISION Le, Argoul 2004
%   Detailed explanation goes here

Df = f/(2*Q);

n = 2*Q^2 - 1/2;
mu_psi = 1/2 * sqrt(1 + 2/(2*n-1));
Dt = mu_psi/Df;
end

