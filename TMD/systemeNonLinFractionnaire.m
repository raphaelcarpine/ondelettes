function [t, X, V, A] = systemeNonLinFractionnaire(variables, equations, variablesNonInertielles, equationsNonInertielles,...
    parametres, valeurParametres, X0, V0, X0NonInertiel, T, nT, dt_integration)
%SYSTEMEQELCONQUE Summary of this function goes here
%   Detailed explanation goes here

vars = variables;
eqs = equations;
ddl1 = length(vars);

param = parametres;
valparam = valeurParametres;
nparam = length(param);
paramStruct = struct;
for kparam = 1:nparam
    paramStruct.(param{kparam}) = valparam(kparam);
end

varsNonI = variablesNonInertielles; % variables avec equa diff 1er ordre seulement
eqsNonI = equationsNonInertielles;
X0nonI = X0NonInertiel;
ddl2 = length(varsNonI);



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


    function str = fracDevRep(str, var, func) % remplace les derives fractionnaires d(alpha)x par func(alpha)
        n = length(str);
        k = length(var);
        position = strfind(str, var);
        hashold = char(35*ones(1, k));
        alphas = {};
        for pos = position
            if pos+k<=n
                l = str(pos+k);
                if l>='a' && l<='z' || l>='A' && l<='Z' || l>='0' && l<='9'
                    continue;
                end
            end
            if pos == 1 || str(pos-1) ~= ')'
                continue;
            end
            
            str(pos:pos+k-1) = hashold;
            
            alpha = '';
            posalpha = pos-2;
            nparenthese = 1;
            while posalpha>1 && nparenthese>0
                alpha = [str(posalpha) alpha];
                if str(posalpha-1) == ')'
                    nparenthese = nparenthese+1;
                elseif str(posalpha-1) == '('
                    nparenthese = nparenthese-1;
                end
                posalpha = posalpha-1;
            end
            
            alphas{end+1} = alpha;
        end
        
        for kalpha = 1:length(alphas)
            alpha = alphas{kalpha};
            str = strrep(str, ['d(' alpha ')' hashold], func(alpha));
        end
    end


    function str = getFracDeriv(alpha, kvar)
        str = ['dt^(-' alpha ') * sum( Aalpha .* Y(' num2str(kvar) ',kt:-1:1))'];
    end



eqs2 = [eqs eqsNonI];

for keqs = 1:length(eqs2)
    for kvar = 1:ddl1
        eqs2{keqs} = varNameRep(eqs2{keqs}, ['d' vars{kvar}], ['Y(' num2str(kvar+ddl1) ')']);
        eqs2{keqs} = varNameRep(eqs2{keqs}, vars{kvar}, ['Y(' num2str(kvar) ')']);
        eqs2{keqs} = fracDevRep(eqs2{keqs}, vars{kvar}, @(alpha) getFracDeriv(alpha, kvar));
    end
    for kvar2 = 1:ddl2
        eqs2{keqs} = varNameRep(eqs2{keqs}, varsNonI{kvar2}, ['Y(' num2str(2*ddl1+kvar2) ')']);
        eqs2{keqs} = fracDevRep(eqs2{keqs}, varsNonI{kvar2}, @(alpha) getFracDeriv(alpha, kvar2));
    end
    for kparam = 1:nparam
        eqs2{keqs} = varNameRep(eqs2{keqs}, param{kparam}, ['paramStruct.' param{kparam}]);
    end
end


Dstring = 'D = @() ['; %string qui va permettre de construire la fonction
for keqs = 1:ddl1
    Dstring = [Dstring 'Y(' num2str(ddl1+keqs) '); '];
end
for keqs = 1:length(eqs2)
    Dstring = [Dstring eqs2{keqs} '; '];
end
Dstring = Dstring(1:end-2);
Dstring = [Dstring '];'];


D = @() 1; % eval ne peut pas déclarer de nouvelle variable

eval(Dstring); % on construit la fonction d'intégration


%% integration

nt = round(T/dt_integration)+1;
dt = T/(nt-1);
Y = zeros(2*ddl1+ddl2, nt);


A

for kt = 1:nt-1
    
end



%%

t = linspace(0, T, T*nT);
Y = interp1(tout, Y, t);
dY = zeros(size(Y'));
for it = 1:length(t)
    dY(:,it) = D(t(it), Y(it,:));
end

X = Y(:,1:ddl1)';
V = Y(ddl1+1:2*ddl1,:);
A = dY(ddl1+1:2*ddl1,:);


end

