% ms = TMDmasseressort(100, 1, @(x,v)0);
% tour1 = Structure(1, 1, @(x,v) 0.0*v, {{ms, 1}});
% tour2 = Structure([1 0;0 100], [2 -1;-1 1], @(x,v) [0.0*v(1);0], {});
% 
% tour1.reponseLibreAvecTMD(0,1, 40);
% tour2.reponseLibreAvecTMD([0 0],[1 0], 40);

pendule = TMDpendule(9.81^2, 1, 9.81, @(theta,omega) 50*omega*(abs(theta)<5*pi/180));
[t, x] = pendule.reponseLibre(pi/2-0.0001, 0, 200);
tempsCaracteristique(t, x(:, 1))

tour = Structure(50, 50, @(x,v) 0.0*v, {{pendule, 1}});

% tour.reponseLibreSansTMD(0, 1, 200);
[t, x] = tour.reponseLibreAvecTMD(0, 1, 200);
tempsCaracteristique(t, x(:, 1))