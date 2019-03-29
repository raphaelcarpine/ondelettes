function systemeQuelconque(variables, equations, parametres, valeurParametres, X0, V0, varargin)
%SYSTEMEQELCONQUE Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

defaultVariablesNonInertielles = {};
defaultEquationsNonInertielles = {};
defaultX0NonInertiel = [];
defaultT = 100;

addRequired(p, 'variables');
addRequired(p, 'equations');
addRequired(p, 'parametres');
addRequired(p, 'valeurParametres');
addRequired(p, 'X0');
addRequired(p, 'V0');
addParameter(p, 'variablesNonInertielles', defaultVariablesNonInertielles);
addParameter(p, 'equationsNonInertielles', defaultEquationsNonInertielles);
addParameter(p, 'X0NonInertiel', defaultX0NonInertiel);
addParameter(p, 'T', defaultT);

parse(p,variables, equations, parametres, valeurParametres, X0, V0, varargin{:})

vars = variables;
eqs = equations;
ddl1 = length(vars);

param = parametres;
valparam = valeurParametres;
nparam = length(param);
S = struct;
for kparam = 1:nparam
    S.(param{kparam}) = valparam(kparam);
end

varsNonI = p.Results.variablesNonInertielles; % variables avec equa diff 1er ordre seulement
eqsNonI = p.Results.variablesNonInertielles;
X0nonI = p.Results.X0NonInertiel;
ddl2 = length(varsNonI);


%integration
T = p.Results.T;
nT = T * 200;


%ondelette
Q = 1;
MaxParallelRidges = 1;
fmin = 0.9;
fmax = 1.1;
NbFreq = 100;

%% construction de la fonction d'intégration

    function str2 = varNameRep(str, var, rep)
    end


% passage de variable symbolique à vecteur
eqs2 = eqs;
for keqs = 1:ddl1
    for kvar = 1:ddl1
        eqs2{keqs} = strrep(eqs2{keqs}, ['d' vars{kvar}], ['Y(' num2str(kvar+ddl1) ')']);
        eqs2{keqs} = strrep(eqs2{keqs}, vars{kvar}, ['Y(' num2str(kvar) ')']);
    end
    for kvar2 = 1:ddl2
        eqs2{keqs} = strrep(eqs2{keqs}, varsNonI{kvar2}, ['Y(' num2str(2*ddl1+kvar2) ')']);
    end
    for kparam = 1:nparam
        eqs2{keqs} = strrep(eqs2{keqs}, param{kparam}, ['S.' param{kparam}]);
    end
end
eqsNonI2 = eqsNonI;
for keqs = 1:ddl2
    for kvar = 1:ddl1
        eqsNonI2{keqs} = strrep(eqsNonI2{keqs}, ['d' vars{kvar}], ['Y(' num2str(kvar+ddl1) ')']);
        eqsNonI2{keqs} = strrep(eqsNonI2{keqs}, vars{kvar}, ['Y(' num2str(kvar) ')']);
    end
    for kvar2 = 1:dd2
        eqsNonI2{keqs} = strrep(eqsNonI2{keqs}, varsNonI{kvar2}, ['Y(' num2str(2*ddl1+kvar2) ')']);
    end
    for kparam = 1:nparam
        eqsNonI2{keqs} = strrep(eqsNonI2{keqs}, param{kparam}, ['S.' param{kparam}]);
    end
end



Dstring = 'D = @(t, Y) ['; %string qui va permettre de construire la fonction
for keqs = 1:ddl1
    Dstring = [Dstring 'Y(' num2str(ddl1+keqs) '); '];
end
for keqs = 1:ddl1
    Dstring = [Dstring eqs2{keqs} '; '];
end
for keqs = 1:ddl2
    Dstring = [Dstring eqsNonI2{keqs} '; '];
end
Dstring = Dstring(1:end-2);
Dstring = [Dstring '];'];


D = @(Y) Y; % eval ne peut pas déclarer de nouvelle variable

%% construction du gui
fig = figure;

paramPan = uipanel('Parent',fig, 'Units', 'normalized');
plotPan = uipanel('Parent',fig, 'Units', 'normalized');
plotparamPan = uipanel('Parent',fig, 'Units', 'normalized');

largeur1 = 0.2;
hauteur1 = 0.1;
marge = 0.01;
paramPan.Position = [marge, marge, largeur1, 1-2*marge];
plotPan.Position = [largeur1+2*marge, hauteur1+2*marge, 1-largeur1-3*marge, 1-hauteur1-3*marge];
plotparamPan.Position = [largeur1+2*marge, marge, 1-largeur1-3*marge, hauteur1];





%%

    function update()
        eval(Dstring); % on construit la fonction d'intégration
        t = linspace(0, T, nT);
        [tout, Xout] = differentialEq([X0; V0; X0nonI], D, T, false, 'output', 'raw');
        X = interp1(tout, Xout, t);
    end

update();

end

