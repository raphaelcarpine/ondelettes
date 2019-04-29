function fig = RegressionMenu(varargin)
%WaveletMenu Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

defaultEq = 'a*x+b';
defaultParam = 'a b';
defaultParam0 = '1 1';
defaultFit = 'y';

paramPrecision = 1e-6;
maxIter = 1e6;

addOptional(p, 'Equation', defaultEq);
addOptional(p, 'Param', defaultParam);
addOptional(p, 'Param0', defaultParam0);
addOptional(p, 'Fit', defaultFit);

parse(p, varargin{:});

eq = p.Results.Equation;
param = p.Results.Param;
param0 = p.Results.Param0;
fit = p.Results.Fit;
%%

fig = figure;
fig.Units = 'characters';
fig.Position(3) = 65;
fig.Position(4) = 22;
fig.MenuBar = 'none';



%% bouton regression et panneaux param et sorties

eqPan = uipanel('Parent',fig, 'Units', 'normalized');
linePan = uipanel('Parent',fig, 'Units', 'normalized');
optionsPan = uipanel('Parent',fig, 'Units', 'normalized');
buttonReg = uicontrol('Parent',fig, 'Units', 'normalized','Style','pushbutton',...
    'String', 'regression');

linePan.Position = [0.02 0.84 0.96 0.14];
eqPan.Position = [0.02 0.37 0.96 0.45];
optionsPan.Position = [0.02 0.15 0.96 0.2];
buttonReg.Position = [0.02 0.02 0.96 0.11];



%% ligne � fitter

line = [];

lines = [];
kline = 0;

highlighted = [];
lineWidth = [];


lineSelect = uicontrol('Parent', linePan, 'Units', 'normalized','Style','togglebutton', 'String', 'select line');
prevBut = uicontrol('Parent', linePan, 'Units', 'normalized','Style','pushbutton', 'String', '<');
nextBut = uicontrol('Parent', linePan, 'Units', 'normalized','Style','pushbutton', 'String', '>');

lineSelect.Position = [0.01, 0.02, 0.48, 0.96];
prevBut.Position = [0.6, 0.1, 0.15, 0.8];
nextBut.Position = [0.75, 0.1, 0.15, 0.8];


    function highlightLine()
        try
            set(highlighted, 'LineWidth', lineWidth);
            lineWidth = get(line, 'LineWidth');
            highlighted = line;
            set(highlighted, 'LineWidth', 3*lineWidth);
        catch
        end
    end

    function selectFunction(selecting)
        if selecting
            lines = findobj('Type', 'line');
            kline = 1;
            line = lines(1:min(1,end));
        else
            line = [];
            lines = [];
        end
        highlightLine();
    end

    function nextprev(next)
        if ~isempty(lines)
            kline = mod(kline + next -1, length(lines)) + 1;
            line = lines(kline);
            highlightLine();
        end
    end

lineSelect.Callback = @(~,~) selectFunction(lineSelect.Value);
prevBut.Callback = @(~,~) nextprev(-1);
nextBut.Callback = @(~,~) nextprev(1);

    function closeReg()
        try
            selectFunction(false);
        catch
        end
        delete(fig);
    end

fig.CloseRequestFcn = @(~,~) closeReg();

%% equation

eqStr = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','text', 'String', 'y = ');
eqEdit = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','edit', 'String', eq);
paramStr = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','text', 'String', 'params :');
paramEdit = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','edit', 'String', param);
param0Str = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','text', 'String', 'params0 :');
param0Edit = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','edit', 'String', param0);
fitStr = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','text', 'String', 'fit :');
fitEdit = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','edit', 'String', fit);

nLign = 4;
marge = 0.02;
h = (1-(nLign+1)*marge)/nLign;
H = h+marge;
eqStr.Position = [0.01, 1-H, 0.2, h];
eqEdit.Position = [0.22, 1-H, 0.77, h];
paramStr.Position = [0.01, 1-2*H, 0.2, h];
paramEdit.Position = [0.22, 1-2*H, 0.77, h];
param0Str.Position = [0.01, 1-3*H, 0.2, h];
param0Edit.Position = [0.22, 1-3*H, 0.77, h];
fitStr.Position = [0.01, 1-4*H, 0.2, h];
fitEdit.Position = [0.22, 1-4*H, 0.77, h];


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
%         line = findobj('Type', 'line');
        if isempty(line)
            ok = false;
        else
%             line = line(1);
            X = get(line, 'XData');
            Y = get(line, 'YData');
            
            X = X(~isnan(Y)); % on enl�ve les valeurs inappropri�es
            Y = Y(~isnan(Y));
            X = X(~isnan(X));
            Y = Y(~isnan(X));
            
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


    function computeReg()
        if ~updateXYAxes()
            warning('no line selected');
            return;
        end
        
        Eq = eqEdit.String;
        Param = paramEdit.String;
        Param0 = param0Edit.String;
        Param = strsplit(Param);
        Param0 = strsplit(Param0);
        for ip = 1:length(Param0)
            Param0{ip} = eval(Param0{ip});
        end
        Param0 = [Param0{:}];
        
        set(param0Edit, 'ForegroundColor', [0.5 0.5 0.5]);
        drawnow;
        
        Fstring = Eq;
        for ip = 1:length(Param)
            Fstring = varNameRep(Fstring, Param{ip}, ['P(' num2str(ip) ')']);
        end
        Fstring = strrep(Fstring, '*', '.*');
        Fstring = strrep(Fstring, '/', './');
        Fstring = strrep(Fstring, '^', '.^');
        
        F = @(P) 0;        
        eval(['F = @(P, x) ' Fstring ';']); 
        
        fitFunction = @(y) y;
        eval(['fitFunction = @(y)' get(fitEdit, 'String') ';']);
        
        S = @(P) fitFunction(F(P, X)) - fitFunction(Y);
        
        lb = ones(size(Param0))*(-inf);
        ub = ones(size(Param0))*inf;
        optionsReg = optimoptions(@lsqnonlin, 'MaxIterations', maxIter,...
            'StepTolerance', paramPrecision, 'MaxFunctionEvaluations', inf, 'FunctionTolerance', 0);
        
        try
            Param1 = lsqnonlin(S, Param0, lb, ub, optionsReg);
        catch error
            warning('did not fit');
            warning(error.message);
            set(param0Edit, 'ForegroundColor', [0 0 0]);
            return
        end
        
        set(param0Edit, 'String', num2str(Param1));
        set(param0Edit, 'ForegroundColor', [0 0 0]);
        
        if optBut.plot.Value
            if length(X) < 1000
                Xplot = linspace(X(1), X(end), 1000);
            else
                Xplot = X;
            end
            
            if optBut.onaxes.Value
                plotAxes = ax;
                hold(plotAxes, 'on');
                onAxesPlots = [onAxesPlots, plot(plotAxes, Xplot, F(Param1, Xplot) .* ones(size(Xplot)), 'r--')];
                hold(plotAxes, 'off');
            else
                plotAxes = axes(figure);
                hold(plotAxes, 'on');
                plot(plotAxes, X, Y, '*');
                plot(plotAxes, Xplot, F(Param1, Xplot) .* ones(size(Xplot)));
                hold(plotAxes, 'off');
            end
        end
        
        lineSelect.Value = false;
        selectFunction(lineSelect.Value);
    end



buttonReg.Callback = @(~,~) computeReg();


end