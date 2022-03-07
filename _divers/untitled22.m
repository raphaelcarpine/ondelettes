dt = 0.02;
N = round(100/dt);
t = dt * (0:N-1);
PsiN = 2;
Psi = 2 + 0.*sin(2*pi*0.1*t);
BetaN = 2;
lambda = 0.1;

coeff = 2 + dt^2*BetaN/lambda;
coeffC = dt^2*BetaN*Psi/lambda;

%%

phi0 = 5;
phiNm1 = 2;

decimalPrec = max([-log10(coeff-2), -log10(coeffC)]) + 5;
prec = ceil(decimalPrec / log10(LongInt.Kmax.double()));
Prec = trunc(LongInt(1), -prec);

A = LongInt.zeros(1, N);
B = LongInt.zeros(1, N);
C = LongInt.zeros(1, N);
A(1) = trunc(LongInt(1), -prec);
B(2) = trunc(LongInt(1), -prec);
coeffP = Prec .* coeff;

for kt = 3:N
%     A(kt) = trunc(coeffP .* A(kt-1), prec) - A(kt-2);
    B(kt) = trunc(coeffP .* B(kt-1), prec) - B(kt-2);
    C(kt) = trunc(coeffP .* C(kt-1), prec) - C(kt-2) - coeffC(kt-1) .* Prec;
end

phiA = nan(1, N);
phiB = nan(1, N);
phiC = nan(1, N);
for kt = 1:N
%     phiA(kt) = divideDouble(A(kt).*B(end) - B(kt).*A(end), B(end).trunc(-prec));
    phiB(kt) = divideDouble(B(kt), B(end));
    phiC(kt) = divideDouble(C(kt).*B(end) - B(kt).*C(end), B(end).trunc(-prec));
end

% phi = phiA*phi0 + phiB*phiNm1 + phiC;


figure;
plot(t, phiC);

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



