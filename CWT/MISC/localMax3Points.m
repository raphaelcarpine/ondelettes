function xmax = localMax3Points(x, y)
%localMax3Points Summary of this function goes here
%   Detailed explanation goes here
if size(x, 1) ~= 3 || size(x, 1) ~= 3
    error('');
end
y1x2x3 = y(1,:).*(x(2,:)-x(3,:));
y2x3x1 = y(2,:).*(x(3,:)-x(1,:));
y3x1x2 = y(3,:).*(x(1,:)-x(2,:));

xmax = (y1x2x3.*(x(2,:)+x(3,:)) + y2x3x1.*(x(3,:)+x(1,:)) + y3x1x2.*(x(1,:)+x(2,:))) ./...
    (2*(y1x2x3 + y2x3x1 + y3x1x2));
end
