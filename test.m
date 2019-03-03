% pendule = TMDpendule(9.81^2, 1, 9.81, @(theta,omega) sign(omega));
% pendule.reponseLibre(0, 2, 500);
% 
% mr = TMDmasseressort(1, 1, @(x, v) 0.01*sign(v));
% mr.reponseLibre(0, 2, 500);
% 
% pendule = TMDpendule(9.81^2, 1, 9.81, @(theta,omega) 50*omega*(abs(theta)<0.1));
% tour = Structure(50, 50, @(x,v) 0.0*v, {{pendule, 1}});
% tour.reponseLibreAvecTMD(0, 1, 200);

% l = 1;
% pendule = TMDpendule(l^2, 1, l, @(theta,omega) 0);
% tour = Structure(50, 50, @(x,v) 10*v, {{pendule, 1}});
% tour.reponseForceeAvecTMD(0, 0, 100, @(t) [1 0], true);

% mr = TMDmasseressort(1, 1, @(x, v) 0.15*v);
% % mr.reponseLibre(0, 1, 200);
% tour = Structure(100, 100, @(x,v) 10*v, {{mr, 1}});
% tour = Structure(100, 100, @(x,v) 0*v, {{mr, 1}});
% tour.reponseLibre(1, 0, 150, true);
% % tour.diagrammeBode(1, 1, 1/(2*pi)*exp(linspace(-1, 1, 1000)), 1000, true);

mr = TMDmasseressort(1, 1, @(x, v) 0.1*v*(abs(x)<0.1));
mr = TMDmasseressort(1, 1, @(x, v) 0.1*sign(v)*(abs(x)<0.1));
mr = TMDmasseressort(1, 1, @(x, v) -0.01*v);
mr = TMDmasseressort(1, 1, @(x, v) 0.05*sign(v)*abs(v)^3);
[t, X] = mr.reponseLibre(0, 1, 200);
x = X(:, 1);
CWT(t, x, 1/(2*pi)*exp(linspace(-1, 1, 100)), 1000);