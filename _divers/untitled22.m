dt = 0.01;
N = round(100/dt);
t = dt * (0:N-1);
PsiN = 2;
Psi = 2 + 0.5*sin(2*pi*0.1*t);
BetaN = 2;
lambda = 1;

coeff = 2 + dt^2*BetaN/lambda;
coeffC = dt^2*BetaN*Psi/lambda;

%%

phi0 = 5;
phiNm1 = 2;

tens = ceil(max([-log10(coeff-2), -log10(coeffC)])) + 10;

A = zeros(1, N);
A(1) = 10^tens;
B = zeros(1, N);
B(2) = 10^tens;
C = zeros(1, N);
coeff = 10^tens * coeff;
coeffC = 10^tens * coeffC;
A = vpi(A);
B = vpi(B);
C = vpi(C);
coeff = vpi(round(coeff));
coeffC = vpi(round(coeffC));
for kt = 3:N
    A(kt) = coeff*A(kt-1)/10^tens - A(kt-2);
    B(kt) = coeff*B(kt-1)/10^tens - B(kt-2);
    C(kt) = coeff*C(kt-1)/10^tens - C(kt-2) - coeffC(kt-1);
end

phi = ((A*B(end) - B*A(end))*phi0 + 10^tens*B*phiNm1 + C*B(end) - B*C(end))/(B(end));
phi = double(phi)/10^tens;


figure;
plot(t, phi);

%%

% A = A/10^tens;
% B = B/10^tens;
% C = C/10^tens;
% 
% M = [dt*BetaN*sum(A.^2) + lambda/dt*sum(dA.^2), dt*BetaN*sum(A.*B) + lambda/dt*sum(dA.*dB);
%     dt*BetaN*sum(A.*B) + lambda/dt*sum(dA.*dB), dt*BetaN*sum(B.^2) + lambda/dt*sum(dB.^2)];
% Y = - [dt*BetaN*sum(A.*C) + lambda/dt*sum(dA.*dC);
%     dt*BetaN*sum(B.*C) + lambda/dt*sum(dB.*dC)];
% 
% phi01 = M\Y;
% 
% % psi = PsiN*[1; 1];
% 
% figure;
% plot(t, [A.', B.', C.'] * [phi01; 1]);



