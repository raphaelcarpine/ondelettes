function Ax = findAxesAudio(ax)
%FINDAXESAUDIO Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0 % test
    ax = gca;
end

%% dialog box

fig = dialog('Name', 'Axis selection menu', 'Units', 'characters', 'WindowStyle', 'normal');
fig.Position(3) = 40;
fig.Position(4) = 7.5;

uicontrol(fig, 'Style', 'text', 'String', 'Axis selection', 'FontWeight', 'bold',...
    'Units', 'characters', 'Position', [0 5.5 40 1]);
txtAxes = uicontrol(fig, 'Style', 'text', 'Units', 'characters', 'Position', [0 3.5 40 1], 'String',...
    ['Click on axis (Fig', num2str(ax.Parent.Number), ' axis selected)']);
okButton = uicontrol(fig, 'Style', 'pushbutton', 'String', 'OK', 'Units', 'characters', 'Position', [7 0.5 12 2]);
cancelButton = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Cancel', 'Units', 'characters', 'Position', [21 0.5 12 2]);



%% axes selection

allAxes = findall(0, 'type', 'axes');
allAxesButtonDownFcns = get(allAxes, 'ButtonDownFcn');
if length(allAxes) <= 1
    allAxesButtonDownFcns = {allAxesButtonDownFcns};
end

set(allAxes, 'ButtonDownFcn', @changeAxis);

Ax = gobjects(1, 0);

    function changeAxis(ax2, ~)
        ax = ax2;
        txtAxes.String = ['Click on axis (Fig', num2str(ax.Parent.Number), ' axis selected)'];
    end


%% callback

cancelButton.Callback = @(~,~) delete(fig);
okButton.Callback = @okCallback;

    function okCallback(~,~)
        Ax = ax;
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

end

