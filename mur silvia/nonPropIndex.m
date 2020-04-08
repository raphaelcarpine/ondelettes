function I = nonPropIndex(shape)
%NONPROPINDEX Summary of this function goes here
%   Detailed explanation goes here

dims = size(shape);
if dims(2) > dims(1)
    shape = transpose(shape);
end

I = norm(imag(shape)) / norm(shape); % sin

% S. Adhikari, Optimal complex modes and an index of damping non-proportionality
shapeOpt = shape * (shape'*real(shape)) / (shape'*shape);
I = norm(shapeOpt - real(shape)) / norm(real(shape));

end

