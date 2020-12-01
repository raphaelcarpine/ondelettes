function dispAverageMenu(AverageNoise, fAverage, multiSignalMode, typeOfAverage, ignoreXLim, SignalUnit)
%DISPAVERAGEMENU Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
    fAverage = 7.6;
    AverageNoise = [pi, 0.2*pi];
    multiSignalMode = false;
    AverageNoise = pi;
    multiSignalMode = true;
    typeOfAverage = 'lin';
    ignoreXLim = false;
    SignalUnit = 'm/s²';
end


%%
TitleMsg = 'Average CWT over time';

message = {};
if multiSignalMode
    switch typeOfAverage
        case 'lin'
            message{end+1} = 'Square root of average of |CWT|²';
        case 'log'
            message{end+1} = 'Square root of logarithmic average of |CWT|²';
    end
    message{end} = [message{end}, sprintf(' for f = %.2f Hz', fAverage)];
    if ignoreXLim
        message{end} = [message{end}, ' (whole signal)'];
    end
    
    message{end+1} = [sprintf(' - all channels: %6e', AverageNoise), ' ', SignalUnit];
else
    switch typeOfAverage
        case 'lin'
            message{end+1} = 'Average of |CWT|';
        case 'log'
            message{end+1} = 'Logarithmic average of |CWT|';
    end
    message{end} = [message{end}, sprintf(' for f = %.2f', fAverage)];
    if ignoreXLim
        message{end} = [message{end}, ' (whole signal)'];
    end
    
    for k_ch = 1:length(AverageNoise)
        message{end+1} = [sprintf(' - channel %d: %6e', [k_ch, AverageNoise(k_ch)]), ' ', SignalUnit];
    end
end

% msgbox(message, TitleMsg, 'help');

%%

fig = figure('Name', TitleMsg, 'numbertitle', 'off');
fig.MenuBar = 'none';
fig.Resize = 'off';
fig.Units = 'characters';
fig.Position(3) = max(60, length(message{1}) + 5);
H = 1.5*length(message) + 3.8;
fig.Position(4) = H;

uicontrol(fig, 'Style', 'text', 'Units', 'character', 'HorizontalAlignment', 'left',...
    'Position', [2, H-2.5, 200, 1.5], 'String', message{1});

for k_ch = 1:length(message)-1
    uicontrol(fig, 'Style', 'text', 'Units', 'character', 'HorizontalAlignment', 'left',...
        'Position', [2, H - 2.5 - 1.5*k_ch, 35, 1.5], 'String', message{k_ch+1});
    copyBut = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'character',...
        'Position', [37, H - 2.5 - 1.5*k_ch, 10, 1.6], 'String', 'copy');
    copyBut.Callback = @(~,~) clipboard('copy', sprintf('%.10e', AverageNoise(k_ch)));
end

okBut = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'character',...
        'Position', [fig.Position(3)/2-5, 0.5, 10, 1.7], 'String', 'OK');
okBut.Callback = @(~,~) delete(fig);
set(okBut, 'Value', true);

end

