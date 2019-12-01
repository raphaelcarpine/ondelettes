function [X, Y] = PlotExtract()
%WaveletMenu Summary of this function goes here
%   Detailed explanation goes here

fig = figure;
fig.Units = 'characters';
fig.Position(3) = 40;
fig.Position(4) = 10;
fig.MenuBar = 'none';



%% bouton regression et panneaux param et sorties

linePan = uipanel('Parent',fig, 'Units', 'normalized');
buttonGet = uicontrol('Parent',fig, 'Units', 'normalized','Style','pushbutton',...
    'String', 'get');

linePan.Position = [0.02 0.51 0.96 0.47];
buttonGet.Position = [0.02 0.02 0.96 0.47];



%% ligne à fitter

line = [];

lines = [];
kline = 0;

highlighted = [];
lineWidth = [];


lineSelect = uicontrol('Parent', linePan, 'Units', 'normalized','Style','togglebutton', 'String', 'select line');
prevBut = uicontrol('Parent', linePan, 'Units', 'normalized','Style','pushbutton', 'String', '<');
nextBut = uicontrol('Parent', linePan, 'Units', 'normalized','Style','pushbutton', 'String', '>');

lineSelect.Position = [0.01, 0.02, 0.48, 0.96];
prevBut.Position = [0.51, 0.01, 0.23, 0.96];
nextBut.Position = [0.76, 0.01, 0.23, 0.96];


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
            
            X = X(~isnan(Y)); % on enlève les valeurs inappropriées
            Y = Y(~isnan(Y));
            X = X(~isnan(X));
            Y = Y(~isnan(X));
            
            ax = get(line, 'Parent');
            ok = true;
        end
    end


%%  

    function getLine()
        if ~updateXYAxes()
            warning('no line selected');
            return;
        end
        closeReg();
    end



buttonGet.Callback = @(~,~) getLine();

waitfor(fig);


end