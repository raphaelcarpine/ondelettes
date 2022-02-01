function signalChannels = getSignalChannels(signalChannels, nbPlots, signalChannelsNames)
%GETSIGNALCHANNELS Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0 % test
    signalChannels = [2 4 5];
    nbPlots = 5;
end


%% display

fig = figure('Name', 'Set signal channels', 'numbertitle', 'off');
fig.Units = 'characters';
fig.Position(3) = 30;
fig.Position(4) = 5 + 1.5*nbPlots;
fig.MenuBar = 'none';
fig.Resize = false;

uicontrol(fig, 'Style', 'text', 'String', 'Set signal channels', 'Units', 'characters',...
    'Position', [2, 3 + 1.5*nbPlots, 30, 1.5], 'HorizontalAlignment', 'left');

channelInput = [];
for k = 1:nbPlots
    channelInput(k) = uicontrol(fig, 'Style', 'checkbox', 'Value', ismember(k, signalChannels), 'String',...
    signalChannelsNames{k}, 'Units', 'characters', 'Position', [3, 3 + 1.5*(nbPlots-k), 30, 1.5],...
    'HorizontalAlignment', 'left');
end

okInput = uicontrol(fig, 'Style', 'pushbutton', 'String','OK', 'Units', 'characters',...
    'Position', [(fig.Position(3)-14)/2, 0.5, 14, 2]);

%% return

okInput.Callback = @okReturn;

    function okReturn(~,~)
        signalChannels = true(1, nbPlots);
        for k_ch = 1:nbPlots
            signalChannels(k_ch) = get(channelInput(k_ch), 'Value');
        end
        signalChannels = find(signalChannels);
        
        delete(fig);
    end

waitfor(fig);
    
end

