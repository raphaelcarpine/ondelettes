function [xmax, zmax] = localMax3Points(x, z)
%localMax3Points Summary of this function goes here
%   Detailed explanation goes here
if size(x, 1) ~= 3 || ~isequal(size(x), size(z))
    warning('size problem');
end

phaseContinuity = false;

y = abs(z);

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
y1 = y(1,:);
y2 = y(2,:);
y3 = y(3,:);
z1 = z(1,:);
z2 = z(2,:);
z3 = z(3,:);

a1 = y1.*(x2-x3);
a2 = y2.*(x3-x1);
a3 = y3.*(x1-x2);

b1 = z1.*(x2-x3);
b2 = z2.*(x3-x1);
b3 = z3.*(x1-x2);



xmax = (a1.*(x2+x3) + a2.*(x3+x1) + a3.*(x1+x2)) ./ (2*(a1+a2+a3));

if phaseContinuity
    ymax = - 1./((x2-x3).*(x3-x1).*(x1-x2)) .*...
        (a1.*(xmax-x2).*(xmax-x3) + a2.*(xmax-x3).*(xmax-x1) + a3.*(xmax-x1).*(xmax-x2));
    
    phimax = nan(size(xmax));
    for k = 1:length(xmax)
        if xmax(k) > x2(k)
            phimax = (xmax(k)-x2(k))/(x3(k)-x2(k)) * angle(z3(k)/z2(k));
        else
            phimax = (xmax(k)-x2(k))/(x1(k)-x2(k)) * angle(z1(k)/z2(k));
        end
    end
    phimax = phimax + angle(z2);
    
    zmax = ymax .* exp(1i*phimax);
else
    zmax = - 1./((x2-x3).*(x3-x1).*(x1-x2)) .*...
        (b1.*(xmax-x2).*(xmax-x3) + b2.*(xmax-x3).*(xmax-x1) + b3.*(xmax-x1).*(xmax-x2));
end

end
