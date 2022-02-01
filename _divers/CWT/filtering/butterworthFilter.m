function Y = butterworthFilter(t, X, f, type, order)
%UNTITLED Filtering of X along dimension 2
%   size(t) == [1, n]
%   size(X) == [d, n];

if ~ismember(type, {'high', 'low'})
    error('unknown type of filter');
end

%%

Nt = length(t);
dt = (t(end)-t(1))/(Nt-1);
fs = 1/dt;

if max( abs( diff(t)/dt -1)) > 1e-4
    error('sampling problem');
end

%%
[B, A] = butter(order, f/(fs/2), type);

Y = filter(B, A, X, [], 2);

end

