function [bode0, bode1] = bodeSystemeLineaire(mu, omega0, omega1, zeta0, zeta1, freqs)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

lambda0 = omega0*zeta0;
lambda1 = omega1*zeta1;

n3 = 2*lambda0;
n2 = omega0^2 + 4*lambda0*lambda1;
n1 = 2*lambda1*omega0^2 + 2*lambda0*omega1^2;
n0 = omega0^2*omega1^2;
d4 = 1;
d3 = 2*lambda0 + 2*lambda1 + 2*lambda1*mu;
d2 = omega0^2 + omega1^2 + 4*lambda0*lambda1 + mu*omega1^2;
d1 = 2*lambda0*omega1^2 + 2*lambda1*omega0^2;
d0 = omega0^2*omega1^2;

H0 = @(p) (n3*p.^3 + n2*p.^2 + n1*p + n0)./(d4*p.^4 + d3*p.^3 + d2*p.^2 + d1*p + d0);


n2 = 4*lambda0*lambda1;
n1 = 2*lambda0*omega1^2 + 2*lambda1*omega0^2;
n0 = omega0^2*omega1^2;
d4 = 1;
d3 = 2*lambda0 + 2*lambda1 + 2*lambda1*mu;
d2 = omega0^2 + omega1^2 + 4*lambda0*lambda1 + mu*omega1^2;
d1 = 2*lambda0*omega1^2 + 2*lambda1*omega0^2;
d0 = omega0^2*omega1^2;

H1 = @(p) (n2*p.^2 + n1*p + n0)./(d4*p.^4 + d3*p.^3 + d2*p.^2 + d1*p + d0);


bode0 = H0(2*1i*pi*freqs);
bode1 = H1(2*1i*pi*freqs);

end


% syms p mu omega0 omega1 lambda0 lambda1;
% 
% M = diag([1, mu]);
% C = 2*mu*[lambda0/mu + lambda1, -lambda1; -lambda1, lambda1];
% K = mu*[omega0^2/mu + omega1^2, -omega1^2; -omega1^2, omega1^2];
% 
% X = (p^2*M+p*C+K)\[p*2*lambda0+omega0^2; 0];
% X(1)
% 
% (
% 2*lambda0*p^3
% + omega0^2*p^2 + 4*lambda0*lambda1*p^2
% + 2*lambda1*omega0^2*p + 2*lambda0*omega1^2*p
% + omega0^2*omega1^2
% )/(
% p^4
% + 2*lambda0*p^3 + 2*lambda1*p^3 + 2*lambda1*mu*p^3
% + omega0^2*p^2 + omega1^2*p^2 + 4*lambda0*lambda1*p^2 + mu*omega1^2*p^2
% + 2*lambda0*omega1^2*p + 2*lambda1*omega0^2*p
% + omega0^2*omega1^2
% )



% syms p mu omega0 omega1 lambda0 lambda1;
% 
% M = diag([1, mu]);
% C = 2*mu*[lambda0/mu + lambda1, -lambda1; -lambda1, lambda1];
% K = mu*[omega0^2/mu + omega1^2, -omega1^2; -omega1^2, omega1^2];
% 
% X = (p^2*M+p*C+K)\[p*2*lambda0+omega0^2; 0];
% X(2)
% 
% 
% (
% 4*lambda0*lambda1*p^2
% + 2*lambda0*omega1^2*p + 2*lambda1*omega0^2*p
% + omega0^2*omega1^2
% )/(
% p^4
% + 2*lambda0*p^3 + 2*lambda1*p^3 + 2*lambda1*mu*p^3
% + omega0^2*p^2 + omega1^2*p^2 + 4*lambda0*lambda1*p^2 + mu*omega1^2*p^2
% + 2*lambda0*omega1^2*p + 2*lambda1*omega0^2*p
% + omega0^2*omega1^2
% )
 















