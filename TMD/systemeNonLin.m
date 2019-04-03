function [t, X, V, A] = systemeNonLin(variables, equations, variablesNonInertielles, equationsNonInertielles,...
    parametres, valeurParametres, X0, V0, X0NonInertiel, T, nT)
%SYSTEMEQELCONQUE Summary of this function goes here
%   Detailed explanation goes here

vars = variables;
eqs = equations; % d2x_i/dt2 = ...
ddl1 = length(vars);

param = parametres;
valparam = valeurParametres;
nparam = length(param);
paramStruct = struct;
for kparam = 1:nparam
    paramStruct.(param{kparam}) = valparam{kparam};
end

varsNonI = variablesNonInertielles; % variables avec equa diff 1er ordre seulement
eqsNonI = equationsNonInertielles; % dx_i/dt = ...
X0nonI = X0NonInertiel;
ddl2 = length(varsNonI);

%%
%integration
solver = @ode45;



%% construction de la fonction d'intégration

    function str = varNameRep(str, old, new) % fonction de remplacement de nom de variable
        n = length(str);
        k = length(old);
        position = strfind(str, old);
        hashold = char(35*ones(1, k));
        for pos = position
            if pos>1
                l = str(pos-1);
                if l>='a' && l<='z' || l>='A' && l<='Z' || l>='0' && l<='9'
                    continue;
                end
            end
            if pos+k<=n
                l = str(pos+k);
                if l>='a' && l<='z' || l>='A' && l<='Z' || l>='0' && l<='9'
                    continue;
                end
            end
            str(pos:pos+k-1) = hashold;
        end
        str = strrep(str, hashold, new);
    end




eqs2 = [eqs eqsNonI];

for keqs = 1:length(eqs2)
    for kvar = 1:ddl1
        eqs2{keqs} = varNameRep(eqs2{keqs}, ['d' vars{kvar}], ['Y(' num2str(kvar+ddl1) ')']);
        eqs2{keqs} = varNameRep(eqs2{keqs}, vars{kvar}, ['Y(' num2str(kvar) ')']);
    end
    for kvar2 = 1:ddl2
        eqs2{keqs} = varNameRep(eqs2{keqs}, varsNonI{kvar2}, ['Y(' num2str(2*ddl1+kvar2) ')']);
    end
    for kparam = 1:nparam
        eqs2{keqs} = varNameRep(eqs2{keqs}, param{kparam}, ['paramStruct.' param{kparam}]);
    end
end


Dstring = 'D = @(t, Y) ['; %string qui va permettre de construire la fonction
for keqs = 1:ddl1
    Dstring = [Dstring 'Y(' num2str(ddl1+keqs) '); '];
end
for keqs = 1:length(eqs2)
    Dstring = [Dstring eqs2{keqs} '; '];
end
Dstring = Dstring(1:end-2);
Dstring = [Dstring '];'];


D = @(Y) Y; % eval ne peut pas déclarer de nouvelle variable

eval(Dstring); % on construit la fonction d'intégration



%% integration

[tout, Xout] = solver(D, [0 T], [X0; V0; X0nonI],...
    odeset('RelTol', 1e-10, 'Stats', 'off', 'MaxStep', 1/nT));


t = linspace(0, T, T*nT);
Y = interp1(tout, Xout, t);
dY = zeros(size(Y'));
for it = 1:length(t)
    dY(:,it) = D(t(it), Y(it,:));
end

X = Y(:,1:ddl1)';
V = Y(:,ddl1+1:2*ddl1)';
A = dY(ddl1+1:2*ddl1,:);       


end

