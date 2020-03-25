function meanX = wmean(x, w)
%WMEAN Weighted mean, removes naN values from X
%   x : vector
%   w : weights

if isequal(size(x), size(w))
    dim = 'all';
elseif size(w, 1) == 1 && size(w, 2) == size(x, 2)
    dim = 1;
elseif size(w, 2) == 1 && size(w, 1) == size(x, 1)
    dim = 2;
else
    error('dimensions');
end

if isequal(dim, 'all')
    w = w(~isnan(x));
    x = x(~isnan(x));
    meanX = sum(x.*w, 'all') / sum(w, 'all');
elseif dim == 1
    meanX = nan(size(x, 1), 1);
    for k = 1:size(x, 1)
        X = x(k, :);
        W = w(~isnan(X));
        X = X(~isnan(X));
        meanX(k, 1) = sum(X.*W) / sum(W);
    end
elseif dim == 2
    meanX = nan(1, size(x, 2));
    for k = 1:size(x, 2)
        X = x(:, k);
        W = w(~isnan(X));
        X = X(~isnan(X));
        meanX(1, k) = sum(X.*W) / sum(W);
    end
end

end

