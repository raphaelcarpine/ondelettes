function poles = polesSystFrac(mu, omega0, omega1, zeta0, zeta1, alpha, precis)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin == 6
    precis = false;
end



maxQ = 100; % denominateur max

%% détermination de p,q tq alpha ~ p/q avec q<=maxQ

p = 1;
q = 1;
for k = 1:maxQ
    l = round(alpha*k);
    if abs(alpha-l/k) < abs(alpha-p/q)
        p = l;
        q = k;
    end
end

%%

lambda0 = zeta0*omega0;
lambda1 = zeta1*omega1^(2-alpha);


Det = @(z)...
1 * z.^4 + ...
2*lambda0 * z.^3 + ...
(omega0^2 + (1+mu)*omega1^2) * z.^2 + ...
2*lambda0*omega1^2 * z + ...
omega0^2*omega1^2 + ...
2*lambda1*(1+mu) * z.^(alpha+2) + ...
4*lambda0*lambda1 * z.^(alpha+1) + ...
2*lambda1*omega0^2 * z.^alpha;


n = max(4*q, 2*p) + 1;
DetPol = zeros(1, n);

DetPol(n-4*q) = 1;
DetPol(n-3*q) = 2*lambda0;
DetPol(n-2*q) = omega0^2 + (1+mu)*omega1^2;
DetPol(n-1*q) = 2*lambda0*omega1^2;
DetPol(n-0*q) = omega0^2*omega1^2;
DetPol(n-2*q-p) = DetPol(n-2*q-p) + 2*lambda1*(1+mu);
DetPol(n-1*q-p) = DetPol(n-1*q-p) + 4*lambda0*lambda1;
DetPol(n-0*q-p) = DetPol(n-0*q-p) + 2*lambda1*omega0^2;

r = roots(DetPol);
R = r.^q;

sortMat = [abs(Det(R)), R];
sortMat = sortrows(sortMat);

R = sortMat(:,2);

if precis
    for k = 1:4
        R(k) = fminsearch(@(z) abs(Det(z)), R(k));
    end
end

poles = R(1:4);

end





% syms p mu omega0 omega1 lambda0 lambda1 alpha;
% 
% M = diag([1, mu]);
% C = 2*[lambda0, 0; 0, 0];
% Calpha = 2*mu*lambda1*[1, -1; -1, 1];
% K = mu*[omega0^2/mu + omega1^2, -omega1^2; -omega1^2, omega1^2];
% 
% det(p^2*M + p*C + p^alpha*Calpha + K)



 

% 1 * p^4
% + 2*lambda0 * p^3
% + (omega0^2 + (1+mu)*omega1^2) * p^2
% + 2*lambda0*omega1^2 * p
% + 2*lambda1*(1+mu) * p^(alpha+2)
% + 4*lambda0*lambda1 * p^(alpha+1)
% + 2*lambda1*omega0^2 * p^alpha
% + omega0^2*omega1^2
 



















