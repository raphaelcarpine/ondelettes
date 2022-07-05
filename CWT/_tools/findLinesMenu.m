function Lines = findLinesMenu(axOrLines, figTitle)
%FINDAXESAUDIO Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0 % test
    axOrLines = gca;
    figTitle = 'test';
end

%% dialog box

fig = dialog('Name', figTitle, 'Units', 'characters', 'WindowStyle', 'normal');
L = 54;
fig.Position(3) = L;
fig.Position(4) = 7.5;

% WinOnTop(fig);

uicontrol(fig, 'Style', 'text', 'String', 'Axis or lines selection', 'FontWeight', 'bold',...
    'Units', 'characters', 'Position', [0 5.5 L 1]);
txtAxes = uicontrol(fig, 'Style', 'text', 'Units', 'characters', 'Position', [0 3.5 L 1], 'String', axOrLinesString());
okButton = uicontrol(fig, 'Style', 'pushbutton', 'String', 'OK', 'Units', 'characters', 'Position', [L/2-13 0.5 12 2]);
cancelButton = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Cancel', 'Units', 'characters', 'Position', [L/2+1 0.5 12 2]);



%% axes selection

allAxes = findall(0, 'type', 'axes');
allAxesButtonDownFcns = get(allAxes, 'ButtonDownFcn');
if length(allAxes) <= 1
    allAxesButtonDownFcns = {allAxesButtonDownFcns};
end
allLines = findall(0, 'type', 'line');
allLinesButtonDownFcns = get(allLines, 'ButtonDownFcn');
if length(allLines) <= 1
    allLinesButtonDownFcns = {allLinesButtonDownFcns};
end

set(allAxes, 'ButtonDownFcn', @changeAxisOrLines);
set(allLines, 'ButtonDownFcn', @changeAxisOrLines);

Lines = gobjects(1, 0);

    function changeAxisOrLines(axOrLines2, ~)
        switch class(axOrLines2)
            case 'matlab.graphics.axis.Axes'
                axOrLines = axOrLines2;
            case 'matlab.graphics.chart.primitive.Line'
                switch class(axOrLines)
                    case 'matlab.graphics.axis.Axes'
                        axOrLines = axOrLines2;
                    case 'matlab.graphics.chart.primitive.Line'
                        if ismember(axOrLines2, axOrLines)
                            axOrLines = axOrLines(axOrLines ~= axOrLines2);
                        else
                            axOrLines(end+1) = axOrLines2;
                        end
                end
        end
        txtAxes.String = axOrLinesString();
    end

    function s = axOrLinesString()
        switch class(axOrLines)
            case 'matlab.graphics.axis.Axes'
                s = sprintf('fig%u axis selected', axOrLines.Parent.Number);
            case 'matlab.graphics.chart.primitive.Line'
                if length(axOrLines) == 1
                    s = '1 line selected';
                else
                    s = sprintf('%u lines selected', length(axOrLines));
                end
                
                linesAxes = [axOrLines.Parent];
                linesFigures = [linesAxes.Parent];
                linesFigNb = [linesFigures.Number];
                linesFigNb = sort(unique(linesFigNb));
                if length(linesFigNb) == 1
                    s = [s, ', fig', num2str(linesFigNb)];
                else
                    s = [s, ', figs', num2str(linesFigNb(1))];
                    for k_nb = 2:length(linesFigNb)
                        s = [s, '&', num2str(linesFigNb(k_nb))];
                    end
                end
        end
        s = ['Click on axis or lines (', s, ')'];
    end


%% callback

cancelButton.Callback = @(~,~) delete(fig);
okButton.Callback = @okCallback;

    function okCallback(~,~)
        switch class(axOrLines)
            case 'matlab.graphics.axis.Axes'
                Lines = findobj(axOrLines, 'Type', 'line');
            case 'matlab.graphics.chart.primitive.Line'
                Lines = axOrLines;
        end
        delete(fig);
    end

fig.KeyReleaseFcn = @enterCallback;

    function enterCallback(~,ev)
        switch ev.Key
            case 'return'
                okCallback;
        end
    end


%% end

waitfor(fig);

for k_ax = 1:length(allAxes)
    set(allAxes(k_ax), 'ButtonDownFcn', allAxesButtonDownFcns{k_ax});
end
for k_li = 1:length(allLines)
    set(allLines(k_li), 'ButtonDownFcn', allLinesButtonDownFcns{k_li});
end

end

