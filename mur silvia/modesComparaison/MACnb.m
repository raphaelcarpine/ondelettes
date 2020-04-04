function mac = MACnb(shape1, shape2)
%MACNB Summary of this function goes here
%   Detailed explanation goes here
mac = abs(shape1'*shape2)^2 / ((shape1'*shape1) * (shape2'*shape2));
end

