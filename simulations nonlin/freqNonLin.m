function f = freqNonLin(f0, func, A)
%FREQNONLIN Summary of this function goes here
%   func = g/m

I = integral(@(theta) func(A*sin(theta)).*sin(theta), 0, 2*pi);

f = f0 + 1/((2*pi)^3*f0*A) * I;

end

