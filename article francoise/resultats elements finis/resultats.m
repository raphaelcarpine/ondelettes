% dims = 22*22*0.1 mm
% verre BK7 : \rho = 4830, E = 53e9, \nu = 0.2


freqs = [457.78, 660.89, 765.84, 1164.0, 1164.0];




e = 2; % epaisseur eau, mm
disp('freqs');
disp(sqrt(0.1*4830/(0.1*4830+e*1000)) * freqs');