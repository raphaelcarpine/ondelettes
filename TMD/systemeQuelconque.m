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
paramStruct = struct;
for kparam = 1:nparam
    paramStruct.(param{kparam}) = valparam(kparam);
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

%% construction de la fonction d'int�gration

    function str = varNameRep(str, old, new)
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


% passage de variable symbolique � vecteur
eqs2 = eqs;
for keqs = 1:ddl1
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
eqsNonI2 = eqsNonI;
for keqs = 1:ddl2
    for kvar = 1:ddl1
        eqsNonI2{keqs} = varNameRep(eqsNonI2{keqs}, ['d' vars{kvar}], ['Y(' num2str(kvar+ddl1) ')']);
        eqsNonI2{keqs} = varNameRep(eqsNonI2{keqs}, vars{kvar}, ['Y(' num2str(kvar) ')']);
    end
    for kvar2 = 1:dd2
        eqsNonI2{keqs} = varNameRep(eqsNonI2{keqs}, varsNonI{kvar2}, ['Y(' num2str(2*ddl1+kvar2) ')']);
    end
    for kparam = 1:nparam
        eqsNonI2{keqs} = varNameRep(eqsNonI2{keqs}, param{kparam}, ['paramStruct.' param{kparam}]);
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


D = @(Y) Y; % eval ne peut pas d�clarer de nouvelle variable

%% construction du gui
fig = figure;

paramPan = uibuttongroup('Parent',fig, 'Units', 'normalized');
plotPan = uipanel('Parent',fig, 'Units', 'normalized');
plotParamPan = uipanel('Parent',fig, 'Units', 'normalized');

largeur1 = 0.2;
hauteur1 = 0.1;
marge = 0.01;
paramPan.Position = [marge, marge, largeur1, 1-2*marge];
plotPan.Position = [largeur1+2*marge, hauteur1+2*marge, 1-largeur1-3*marge, 1-hauteur1-3*marge];
plotParamPan.Position = [largeur1+2*marge, marge, 1-largeur1-3*marge, hauteur1];

%% params

paramEdits = struct;
% marge = 0.01;
% for kparam=1:nparam
%     uicontrol('Parent',paramPan, 'Units', 'normalized',...
%         'Style','text', 'String', param{kparam},...
%         'Position', [marge, (nparam-kparam)/nparam+marge, 1/2-3/2*marge, 1/nparam-2*marge]);
%     paramEdits.(param{kparam}) = uicontrol('Parent',paramPan, 'Units', 'normalized',...
%         'Style','edit', 'String', num2str(valparam(kparam)), ...
%         'Position', [1/2-1/2*marge, (nparam-kparam)/nparam+marge, 1/2-3/2*marge, 1/nparam-2*marge]);
% end
marge = 0.05;
paramPan.Units = 'characters';
H = paramPan.Position(4);

for kparam=1:nparam
    paramStr = uicontrol('Parent',paramPan, 'Units', 'character',...
        'Style','text', 'String', param{kparam},...
        'Position', [marge, H-2.5*kparam+marge, 10-3/2*marge, 2-2*marge]);
    paramEdits.(param{kparam}) = uicontrol('Parent',paramPan, 'Units', 'character',...
        'Style','edit', 'String', num2str(valparam{kparam}), ...
        'Position', [marge+10, H-2.5*kparam+marge, 10-3/2*marge, 2-2*marge]);
    paramStr.Units = 'normalized';
    paramEdits.(param{kparam}).Units = 'normalized';
end
paramPan.Units = 'normalized';


for kparam = 1:nparam
    paramEdits.(param{kparam}).Callback = @(~,~) update();
end

%% plotparam

plotParam = {'plot', 'T'};
valplotParam = {'x', num2str(T)};
nplotParam = length(plotParam);

plotParamEdits = struct;
plotParamStrings = struct;
marge = 0.03;
for kplotParam=1:nplotParam
    plotParamStrings.(plotParam{kplotParam}) = uicontrol('Parent',plotParamPan, 'Units', 'normalized',...
        'Style','text', 'String', plotParam{kplotParam});
    plotParamEdits.(plotParam{kplotParam}) = uicontrol('Parent',plotParamPan, 'Units', 'normalized',...
        'Style','edit', 'String', valplotParam{kplotParam});
end

plotParamStrings.plot.Position = [marge, marge, 0.1, 1-2*marge];
plotParamEdits.plot.Position = [0.1+2*marge, marge, 0.1, 1-2*marge];

plotParamStrings.T.Position = [0.8-2*marge, marge, 0.1, 1-2*marge];
plotParamEdits.T.Position = [0.9-marge, marge, 0.1, 1-2*marge];


for kplotParam = 1:nplotParam
    plotParamEdits.(plotParam{kplotParam}).Callback = @(~,~) update();
end

%% plot

for kvar = 1:ddl1
    axesPlots(kvar) = subplot(ddl1, 1, kvar, axes('Parent', plotPan));
    axesPlots(kvar).XGrid = 'on';
    axesPlots(kvar).YGrid = 'on';
    plots(kvar) = plot(nan, 'Parent', axesPlots(kvar));
    ylabel(axesPlots(kvar), vars{kvar});
end

xlabel(axesPlots(ddl1), 't');
% linkaxes(axesPlots,'x');


%%

    function update()
        % on masque l'ancien plot
        for plt=plots
            plt.Color(4) = 0.5;
        end
        drawnow;
        
        % mise � jour des parametres
        for kp=1:nparam
            paramStruct.(param{kp}) = eval(get(paramEdits.(param{kp}), 'String'));
        end
        T = eval(get(plotParamEdits.T, 'String'));
        plotValue = get(plotParamEdits.plot, 'String');
        
        % evaluation de la solution
        eval(Dstring); % on construit la fonction d'int�gration
        t = linspace(0, T, nT);
        [tout, Xout] = differentialEq([X0; V0; X0nonI], D, T, false, 'output', 'raw');
        X = interp1(tout, Xout, t);
        if plotValue == 'x'
            Xplot = X(:,1:ddl1)';
        elseif plotValue == 'v'
            Xplot = X(:,ddl1+1:2*ddl1)';
        elseif plotValue == 'a'
            dX = zeros(size(X'));
            for it = 1:length(t)
                dX(:,it) = D(t(it), X(it,:));
            end
            Xplot = dX(ddl1+1:2*ddl1,:);
        end
        
        %affichage
        for kv = 1:ddl1
            set(plots(kv), 'XData', t);
            set(plots(kv), 'YData', Xplot(kv,:));
            plots(kv).Color(4) = 1;
        end
        %linkaxes(axesPlots,'x');
    end

%%
update();

WaveletMenu('WaveletPlot', plots, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);

end

