function Axes = findLinkedAxesAudio()
%FINDAXESAUDIO Summary of this function goes here
%   Detailed explanation goes here

%% dialog box

fig = dialog('Name', 'Linked axes selection menu', 'Units', 'characters', 'WindowStyle', 'normal');
fig.Position(3) = 40;
fig.Position(4) = 7.5;

uicontrol(fig, 'Style', 'text', 'String', 'Linked axes selection', 'FontWeight', 'bold',...
    'Units', 'characters', 'Position', [0 5.5 40 1]);
txtAxes = uicontrol(fig, 'Style', 'text', 'String', 'Click on axes (0 axes selected)',...
    'Units', 'characters', 'Position', [0 3.5 40 1]);
okButton = uicontrol(fig, 'Style', 'pushbutton', 'String', 'OK', 'Units', 'characters', 'Position', [7 0.5 12 2]);
cancelButton = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Cancel', 'Units', 'characters', 'Position', [21 0.5 12 2]);



%% axes selection

allAxes = findall(0, 'type', 'axes');
allAxesButtonDownFcns = get(allAxes, 'ButtonDownFcn');
if length(allAxes) <= 1
    allAxesButtonDownFcns = {allAxesButtonDownFcns};
end

set(allAxes, 'ButtonDownFcn', @addAxesToList);

Axes = gobjects(1, 0);
Axes0 = gobjects(1, 0);

    function addAxesToList(ax, ~)
        if ~ismember(ax, Axes0)
            Axes0(end+1) = ax;
        else
            Axes0 = Axes0(Axes0 ~= ax);
        end
        
        if length(Axes0) == 1
            txtAxes.String = 'Click on axes (1 axis selected)';
        else
            txtAxes.String = sprintf('Click on axes (%u axes selected)', length(Axes0));
        end
    end


%% callback

cancelButton.Callback = @(~,~) delete(fig);
okButton.Callback = @okCallback;

    function okCallback(~,~)
        Axes = Axes0;
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

