function x = meanSmooth(x, N)
%SMOOTH Summary of this function goes here
%   Detailed explanation goes here

X = x;

for k = 1:N
    X = X + [x(1)*ones(1, k), x(1:end-k)];
    X = X + [x(1+k:end), x(end)*ones(1, k)];
end

x = X/(2*N+1);

end

