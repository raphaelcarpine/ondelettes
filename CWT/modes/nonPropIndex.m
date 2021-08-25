function In = nonPropIndex(shape)
%NONPROPINDEX Summary of this function goes here
%   Detailed explanation goes here

if ~iscolumn(shape)
    shape = shape.';
end

if isreal(shape)
    warning('real shape');
end

if abs(shape.'*shape - 1) > 1e-6
    warning('shape^T shape =\= 1');
    shape = shape / sqrt(shape.'*shape);
end

In = (imag(shape)'*imag(shape)) / (shape'*shape);

end

