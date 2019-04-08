function plots = systemeQuelconque(variables, equations, parametres, valeurParametres, X0, V0, fractionnaire, varargin)
%SYSTEMEQELCONQUE Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

defaultVariablesNonInertielles = {};
defaultEquationsNonInertielles = {};
defaultX0NonInertiel = [];
defaultT = 100;
defaultnT = 100;
defaultdt = 1e-3;

addRequired(p, 'variables');
addRequired(p, 'equations');
addRequired(p, 'parametres');
addRequired(p, 'valeurParametres');
addRequired(p, 'X0');
addRequired(p, 'V0');
addRequired(p, 'lineaire');
addParameter(p, 'variablesNonInertielles', defaultVariablesNonInertielles);
addParameter(p, 'equationsNonInertielles', defaultEquationsNonInertielles);
addParameter(p, 'X0NonInertiel', defaultX0NonInertiel);
addParameter(p, 'T', defaultT);
addParameter(p, 'nT', defaultnT);
addParameter(p, 'dt', defaultdt);

parse(p,variables, equations, parametres, valeurParametres, X0, V0, fractionnaire, varargin{:})

vars = variables;
eqs = equations;
ddl1 = length(vars);

param = parametres;
valparam = valeurParametres;
nparam = length(param);

varsNonI = p.Results.variablesNonInertielles; % variables avec equa diff 1er ordre seulement
eqsNonI = p.Results.equationsNonInertielles;
X0NonI = p.Results.X0NonInertiel;

%%
%integration
T = p.Results.T;
nT = p.Results.nT;

dt_fractionnaire = p.Results.dt;


%% construction du gui
fig = figure;

paramPan = uibuttongroup('Parent',fig, 'Units', 'normalized');
plotPan = uipanel('Parent',fig, 'Units', 'normalized');
plotParamPan = uipanel('Parent',fig, 'Units', 'normalized');

largeur1 = 0.2;
hauteur1 = 0.07;
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
if fractionnaire
    plotParam{end+1} = 'dt';
    valplotParam{end+1} = num2str(dt_fractionnaire);
end
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

if fractionnaire
    plotParamStrings.dt.Position = [0.4-1/2*marge, marge, 0.1, 1-2*marge];
    plotParamEdits.dt.Position = [0.5+1/2*marge, marge, 0.1, 1-2*marge];
end


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
        
        % mise à jour des parametres
        for kp=1:nparam
            valparam{kp} = eval(get(paramEdits.(param{kp}), 'String'));
        end
        T = eval(get(plotParamEdits.T, 'String'));
        plotValue = get(plotParamEdits.plot, 'String');
        
        % evaluation de la solution
        if fractionnaire
            dt_fractionnaire = eval(get(plotParamEdits.dt, 'String'));
            [t, X, V, A] = systemeNonLinFractionnaire(vars, eqs, varsNonI, eqsNonI,...
                param, valparam, X0, V0, X0NonI, T, nT, dt_fractionnaire);
        else
            [t, X, V, A] = systemeNonLin(vars, eqs, varsNonI, eqsNonI,...
                param, valparam, X0, V0, X0NonI, T, nT);
        end
        
        if plotValue == 'x'
            Xplot = X;
        elseif plotValue == 'v'
            Xplot = V;
        elseif plotValue == 'a'
            Xplot = A;
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

end

