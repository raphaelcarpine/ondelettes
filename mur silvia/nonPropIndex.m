function I = nonPropIndex(shape)
%NONPROPINDEX Summary of this function goes here
%   Detailed explanation goes here

dims = size(shape);
if dims(2) > dims(1)
    shape = transpose(shape);
end

I = sqrt(imag(shape)'*imag(shape) / (real(shape)'*real(shape))); % tan

I = sqrt(imag(shape)'*imag(shape) / (shape'*shape) ); % sin

% S. Adhikari, Optimal complex modes and an index of damping non-proportionality
% shapeOpt = 
% I = 0;

end

