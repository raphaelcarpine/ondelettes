function [X, Y, K, stdY] = averagingScatter2(x, y, k)
%AVERAGINGSCATTER Averaging scatter plot, >= k points for average
%   Detailed explanation goes here

if k > length(x)
    error('k > length(x)');
elseif k <= 0
    error('k <= 0');
end

if ~iscolumn(x) && ~isrow(x)
    error('x matrix');
elseif ~iscolumn(y) && ~isrow(y)
    error('y matrix');
elseif length(x) ~= length(y)
    error('length(x) ~= length(y)');
end


%%

[x, I] = sort(x);
y = y(I);


N = floor(length(x)/k);

k = floor(length(x)/N);
r = mod(length(x), k);
randomK = 1:N;
for i = 1:(N-r)
    randomK(randi(length(randomK))) = [];
end
K = k * ones(1, N);
K(randomK) = k+1;

%%

X = nan(1, N);
Y = nan(1, N);
stdY = nan(1, N);

for n = 1:N
    X(n) = mean(x(1:K(n)));
    Y(n) = mean(y(1:K(n)));
    stdY(n) = std(y(1:K(n)));
    
    x = x(K(n)+1:end);
    y = y(K(n)+1:end);
end

end

