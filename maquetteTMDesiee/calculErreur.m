Dm = 2.722;

f0(1) = 1.02539;
f1(1) = 1.00098;
f0(2) = 1.0267;
f1(2) = 1;

T(1) = 4096/50;
T(2) = 600;
DDm = 0.002;


dmdDm = @(k) 1/((f0(k)/f1(k))^2-1);
dmdf0 = @(k) -2*f0(k)/f1(k)^2 * Dm/((f0(k)/f1(k))^2-1)^2;
dmdf1 = @(k) 2*f0(k)^2/f1(k)^3 * Dm/((f0(k)/f1(k))^2-1)^2;

sigmam = @(k) sqrt((dmdDm(k)*DDm/(2*sqrt(3)))^2 +...
    (dmdf0(k)*1/(2*sqrt(3)*T(k)))^2 + (dmdf1(k)*1/(2*sqrt(3)*T(k)))^2);

sigmam(1)
sigmam(2)