function fig = RegressionMenu(varargin)
%WaveletMenu Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

defaultEq = 'a*x+b';
defaultParam = 'a b';
defaultParam0 = '1 1';

paramPrecision = 1e-6;
maxIter = 1e6;

addOptional(p, 'Equation', defaultEq);
addOptional(p, 'Param', defaultParam);
addOptional(p, 'Param0', defaultParam0);

parse(p, varargin{:});

eq = p.Results.Equation;
param = p.Results.Param;
param0 = p.Results.Param0;
%%

fig = figure;
fig.Units = 'characters';
fig.Position(3) = 65;
fig.Position(4) = 15;
fig.MenuBar = 'none';



%% bouton regression et panneaux param et sorties
buttonReg = uicontrol('Parent',fig, 'Units', 'normalized','Style','pushbutton',...
    'String', 'regression');
optionsPan = uipanel('Parent',fig, 'Units', 'normalized');
eqPan = uipanel('Parent',fig, 'Units', 'normalized');

buttonReg.Position = [0.02 0.02 0.96 0.15];
optionsPan.Position = [0.02 0.19 0.96 0.33];
eqPan.Position = [0.02 0.54 0.96 0.42];

%% equation

eqStr = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','text', 'String', 'y = ');
eqEdit = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','edit', 'String', eq);
paramStr = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','text', 'String', 'params :');
paramEdit = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','edit', 'String', param);
param0Str = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','text', 'String', 'params0 :');
param0Edit = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','edit', 'String', param0);

eqStr.Position = [0.01, 0.67, 0.2, 0.32];
eqEdit.Position = [0.22, 0.67, 0.77, 0.32];
paramStr.Position = [0.01, 0.34, 0.2, 0.32];
paramEdit.Position = [0.22, 0.34, 0.77, 0.32];
param0Str.Position = [0.01, 0.01, 0.2, 0.32];
param0Edit.Position = [0.22, 0.01, 0.77, 0.32];


%% options

options = {'plot', 'onaxes'};
nopt = length(options);

optionsStr = struct;
optionsStr.plot = 'afficher';
optionsStr.onaxes = 'sur les axes';

optBut = struct;
for kopt = 1:length(options)
    opt = options{kopt};
    optBut.(opt) = uicontrol('Parent', optionsPan, 'Units', 'normalized','Style','checkbox',...
        'String', optionsStr.(opt), 'Position', [0.01, (nopt-kopt)/nopt, 0.98, 1/nopt]);
end
optBut.plot.Value = true;


onAxesPlots = [];
deleteBut = uicontrol('Parent', optionsPan, 'Units', 'normalized','Style','pushbutton',...
    'String', 'delete plots', 'Position', [0.6, (nopt-2)/nopt, 0.39, 1/nopt]);

    function deletePlots()
        delete(onAxesPlots);
        onAxesPlots = [];
    end

deleteBut.Callback = @(~,~) deletePlots();


%% donnees

X = nan;
Y = nan;
ax = 0;

    function ok = updateXYAxes()
        line = findobj('Type', 'line');
        if isempty(line)
            ok = false;
        else
            line = line(1);
            X = get(line, 'XData');
            Y = get(line, 'YData');
            ax = get(line, 'Parent');
            ok = true;
        end
    end


%% construction de la fonction


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



%%


    function show()
        if ~updateXYAxes()
            warning('no line selected');
            return;
        end
        
        Eq = eqEdit.String;
        Param = paramEdit.String;
        Param0 = param0Edit.String;
        Param = strsplit(Param);
        Param0 = str2double(strsplit(Param0));
        
        Fstring = Eq;
        for ip = 1:length(Param)
            Fstring = varNameRep(Fstring, Param{ip}, ['P(' num2str(ip) ')']);
        end
        Fstring = strrep(Fstring, '*', '.*');
        Fstring = strrep(Fstring, '/', './');
        Fstring = strrep(Fstring, '^', '.^');
        
        F = @(P) 0;        
        eval(['F = @(P, x) ' Fstring ';']); 
        
        S = @(P) F(P, X)-Y;
        
        lb = ones(size(Param0))*(-inf);
        ub = ones(size(Param0))*inf;
        optionsReg = optimoptions(@lsqnonlin, 'MaxIterations', maxIter,...
            'StepTolerance', paramPrecision, 'MaxFunctionEvaluations', inf);
        Param1 = lsqnonlin(S, Param0, lb, ub, optionsReg);
        
        set(param0Edit, 'String', num2str(Param1));
        
        if optBut.plot.Value
            if length(X) < 1000
                Xplot = linspace(X(1), X(end), 1000);
            else
                Xplot = X;
            end
            
            if optBut.onaxes.Value
                plotAxes = ax;
                hold(plotAxes, 'on');
                onAxesPlots = [onAxesPlots, plot(plotAxes, Xplot, F(Param1, Xplot), '--')];
                hold(plotAxes, 'off');
            else
                plotAxes = axes(figure);
                hold(plotAxes, 'on');
                plot(plotAxes, X, Y, '*');
                plot(plotAxes, Xplot, F(Param1, Xplot));
                hold(plotAxes, 'off');
            end
        end
    end



buttonReg.Callback = @(~,~) show();


end