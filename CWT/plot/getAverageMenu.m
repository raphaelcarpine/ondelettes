function [f, ignoreXLim, typeOfAverage, fig] = getAverageMenu(typeOfAverageInit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fInit = [];
ignoreXLimInit = true;

if nargin == 0
    typeOfAverageInit = 'lin'; % 'log'
end

%% input

fig = figure('Name', 'Mean CWT Menu', 'numbertitle', 'off');
fig.Units = 'characters';
fig.Position(3) = 52;
fig.Position(4) = 12.5;
fig.MenuBar = 'none';
fig.Resize = false;


% freq
uicontrol('Parent', fig, 'Style', 'text', 'String', 'Frequency [Hz]:',...
    'Units', 'characters', 'Position', [2, 10, 17, 1.5]);
fInput = uicontrol('Parent', fig, 'Style', 'edit', 'String', fInit,...
    'Units', 'characters', 'Position', [19, 10, 14, 1.7]);

% type of average
bg = uibuttongroup('Parent', fig, 'Units', 'characters', 'Position', [1, 5.5, 50, 4], 'Visible', 'on');
linInput = uicontrol(bg, 'Style', 'radiobutton', 'String', 'Time average',...
    'Units', 'characters', 'Position', [0.5, 2.1, 30, 1.7]);
logInput = uicontrol(bg, 'Style', 'radiobutton', 'String', 'Logarithmic time average',...
    'Units', 'characters', 'Position', [0.5, 0.2, 30, 1.7]);
switch typeOfAverageInit
    case 'lin'
        set(linInput, 'Value', true);
    case 'log'
        set(logInput, 'Value', true);
end

% ignoreXLim
ignoreXLimInput = uicontrol('Parent', fig, 'Style', 'checkbox', 'String', 'Ignore XLim (take whole signal)',...
    'Value', ignoreXLimInit, 'Units', 'characters', 'Position', [2, 3.5, 40, 1.5]);

% ok button
okInput = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'OK',...
    'Units', 'characters', 'Position', [21, 0.5, 10, 2.2]);

% focus
uicontrol(fInput);

%% return

set(okInput, 'Callback', @(~,~) okReturn());

    function okReturn()
        set(okInput, 'Value', true);
        set(okInput, 'Selected', true);
        
        f = str2double(get(fInput, 'String'));
        ignoreXLim = get(ignoreXLimInput, 'Value');
        
        if get(linInput, 'Value')
            typeOfAverage = 'lin';
        else
            typeOfAverage = 'log';
        end
        
        set(fig, 'UserData', 'ok');
        
        set(fig, 'pointer', 'watch');
        drawnow;
    end

% shortcut
    function keyboardShortcut(~, event)
        switch event.Key
            case 'return'
                okReturn();
        end
    end
set(fig, 'KeyReleaseFcn', @keyboardShortcut);
set(fInput, 'KeyReleaseFcn', @keyboardShortcut);
% set(fInput, 'Callback', @(~,~) okReturn);

%%
f = fInit;
ignoreXLim = ignoreXLimInit;
typeOfAverage = typeOfAverageInit;

%%
waitfor(fig, 'UserData');

end

