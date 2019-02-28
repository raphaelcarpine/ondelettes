function T = tempsCaracteristique(t, x, x0)
%tempsCaracteristique Summary of this function goes here
%   Detailed explanation goes here
x = x - mean(x);

switch nargin
    case 2
        x0 = x(1);
        x0 = max(abs(x));
end

T = 0;
for j = 1:length(t)-1
    T = T + (x(j)^2+x(j+1)^2)/2 * (t(j+1)-t(j));
end
T = 4*T/x0^2;
end

