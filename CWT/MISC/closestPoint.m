function [k, r] = closestPoint(L, x)
%CLOSESTPOINT Summary of this function goes here
%   L: sorted list
%   x: point, L(1) <= x <= L(end)
%   L(k) <= x < L(k+1)
%   r = (x - L(k1)) / (L(k2) - L(k1))

throxEr = true;

if x < L(1) && throxEr
    error('x < L(1)');
elseif x > L(end) && throxEr
    error('x > L(end)');
end

%%

k1 = 1;
k2 = length(L);

while k1 < k2-1
    km = round((k1+k2)/2);
    if x < L(km)
        k2 = km;
    else
        k1 = km;
    end
end

%%

k = k1;
r = (x - L(k1)) / (L(k2) - L(k1));


end

