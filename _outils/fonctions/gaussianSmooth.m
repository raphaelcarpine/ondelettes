function x = gaussianSmooth(x, N)
%SMOOTH Summary of this function goes here
%   Detailed explanation goes here

for k = 1:N^2
    x = ([x(1), x(1:end-1)] + x + [x(2:end), x(end)])/3;
end

end

