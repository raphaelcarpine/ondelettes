function freqs = getFreqs3ddl(f1, f2, f3, m1, m2, m3, Deltam)
%GETPOLES3DDL Summary of this function goes here
%   Detailed explanation goes here

w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
Dm = Deltam;


detM = [ - Dm*m1*m2 - Dm*m1*m3 - Dm*m2*m3 - m1*m2*m3, Dm*m1*m2*w1^2 + Dm*m1*m2*w2^2 + Dm*m1*m3*w1^2 + Dm*m1*m3*w3^2 + Dm*m2*m3*w2^2 + Dm*m2*m3*w3^2 + m1*m2*m3*w1^2 + m1*m2*m3*w2^2 + m1*m2*m3*w3^2, - Dm*m1*m2*w1^2*w2^2 - Dm*m1*m3*w1^2*w3^2 - Dm*m2*m3*w2^2*w3^2 - m1*m2*m3*w1^2*w2^2 - m1*m2*m3*w1^2*w3^2 - m1*m2*m3*w2^2*w3^2, m1*m2*m3*w1^2*w2^2*w3^2];


racines = roots(detM);

racines = sort(racines);
racines = transpose(racines);

freqs = sqrt(racines) / (2*pi);

end



% syms x w1 w2 w3 m1 m2 m3 Dm
% 
% M = [m1, 0, 0;
%     0, m2, 0;
%     0, 0, m3];
% 
% Om2 = [w1^2, 0, 0;
%     0, w2^2, 0;
%     0, 0, w3^2];
% 
% J = [1, 1, 1;
%     1, 1, 1;
%     1, 1, 1];
% 
% d = det(-x*M + M*Om2 - x*Dm*J);
% 
% Poly = coeffs(d, x);
% 
% fliplr(Poly)