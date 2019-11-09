f0 = 6.653;
f1 = 6.617;
zeta1 = 0.09;
mu = 0.01472;

w0 = 2*pi*f0;
w1 = 2*pi*f1;
z1 = zeta1;


% calcul formel jacobien réduit

syms p [1 4] % poles : wa, wb, la, lb
syms P [1 4] % poles : Oma^2, Omb^2, la, lb
syms a [1 4] % a : a0, a1, a2, a3
syms s [1 4] % systeme : w0, w1, zeta1, mu



P = [p(1)^2 + p(3)^2, p(2)^2 + p(4)^2, p(3), p(4)];

a = [P(1)*P(2),...
    2 * (P(1)*P(4) + P(2)*P(3)),...
    P(1) + P(2) + 4*P(3)*P(4),...
    2*(P(3) + P(4))];

s = [sqrt(a(3) - a(1)*a(4)/a(2)),...
    sqrt(a(1)*a(2) / (a(2)*a(3) - a(1)*a(4))),...
    a(2) / (2 * (a(3) - a(1)*a(4)/a(2)) * sqrt(a(1)*a(2) / (a(2)*a(3) - a(1)*a(4)))),...
    a(4)/a(2) * (a(3) - a(1)*a(4)/a(2)) - 1];


J = jacobian(s, p);


syms Cs [4 4]
syms Cp [4 4]

for i = 1:4
    for j = 1:4
        Cs(i, j) = 0;
        Cp(i, j) = 0;
    end
end

for k = 1:4
    Cs(k, k) = s(k);
    Cp(k, k) = p(k);
end


Jreduit = Cs\J*Cp;



% calcul numérique

det = [1, 2*(1+mu)*z1*w1, w0^2+(1+mu)*w1^2, 2*z1*w0^2*w1, w0^2*w1^2];

r = roots(det);

paramP = [imag(r(1)), imag(r(3)), -real(r(1)), -real(r(3))];

M = subs(Jreduit, p, paramP);
M = eval(M);

disp(M);


% 0.5803    0.4145   -0.0228    0.0279
% 0.4185    0.5846    0.0240   -0.0271
% -0.6186   -0.3762    0.5419    0.4529
% 13.7991  -14.3676    0.2857    0.2828





