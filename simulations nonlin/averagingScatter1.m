function [X, Y, K, stdY] = averagingScatter1(x, y, dx)
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
Y0 = cell(1, N);

for kx = 1:length(x)
    n = floor(x(kx)/dx) - ni + 1;
    Y0{n} = [Y0{n}, y(kx)];
end

Y = nan(1, N);
K = nan(1, N);
stdY = nan(1, N);

for n = 1:N
    Y(n) = mean(Y0{n});
    K(n) = length(Y0{n});
    stdY(n) = std(Y0{n});
end

end

