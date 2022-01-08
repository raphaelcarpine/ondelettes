function [X, Y, stdY, K, Xlims] = averagingScatter1(x, y, dx)
%AVERAGINGSCATTER Summary of this function goes here
%   Detailed explanation goes here


if ~iscolumn(x) && ~isrow(x)
    error('x matrix');
elseif ~iscolumn(y) && ~isrow(y)
    error('y matrix');
elseif length(x) ~= length(y)
    error('length(x) ~= length(y)');
end

%%

n = floor(x/dx);
ni = min(n);
nf = max(n);
N = nf - ni + 1;

X = ((ni:nf)+1/2) * dx;
Xlims = (ni:nf+1) * dx;
Y = zeros(1, N);
K = zeros(1, N);
stdY = zeros(1, N);

for kx = 1:length(x)
    n = floor(x(kx)/dx) - ni + 1;
    K(n) = K(n) + 1;
    Y(n) = Y(n) + y(kx);
    stdY(n) = stdY(n) + y(kx)^2;
end

Y = Y./K;
stdY = sqrt(stdY./K - Y.^2);

end

