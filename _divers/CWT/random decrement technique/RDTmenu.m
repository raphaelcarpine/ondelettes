function [trdt, Xrdt, axRdt] = RDTmenu(t, X, signalChannels, signalUnit)
%RDTMENU Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
    t = linspace(0, 100, 10000);
    X = [2; -1] * sin(t) + [1; 0];
    signalChannels = 2 + (1:size(X, 1));
end

removeMean0 = true;
lagTime0 = 0;
thresholdValue0 = 1;
refCh0 = 1;

%% input array size

% array size
if size(t, 1) == 1 && size(t, 2) == size(X, 2)
    % ok
elseif size(t, 1) == size(X, 1) && size(t, 2) == size(X, 2)
    dt = mean(diff(t(1, :)));
    for k_t = 2:size(t, 1)
        if max(abs(t(k_t, :) - t(1, :))) > dt * 1e-3
            error('different time arrays');
        end
    end
    t = t(1, :);
else
    error(['array size problem (', num2str(size(t)), ' & ', num2str(size(X)), ')']);
end

% time step
if any(abs(diff(t)/mean(diff(t)) - 1) > 1e-3)
    error('non-constant time step');
end

%% standard dev

StdDev = std(X, 0, 2);


%% figure

fig = figure('Name', 'RDT Menu', 'numbertitle', 'off');
fig.Units = 'characters';
fig.Position(1) = fig.Position(1) + fig.Position(3);
fig.Position(3) = 52;
fig.Position(4) = 16.5;
fig.MenuBar = 'none';
fig.Resize = false;


%% display

plotRDTtimes = PlotRDTtimes(t, X, removeMean0, refCh0, signalChannels, signalUnit);

%% menu

% reference channel
refChNames = cell(1, size(X, 1));
for k_ch = 1:size(X, 1)
    refChNames{k_ch} = sprintf('channel %u', signalChannels(k_ch));
end
uicontrol('Parent', fig, 'Style', 'text', 'String', 'Reference channel:',...
    'Units', 'characters', 'Position', [2, 14.2, 22, 1.5], 'HorizontalAlignment', 'left');
refChInput = uicontrol('Parent', fig, 'Style', 'popupmenu', 'String', refChNames,...
    'Units', 'characters', 'Position', [23, 14.2, 15, 1.7], 'Value', refCh0);

% remove mean
removeMeanInput = uicontrol('Parent', fig, 'Style', 'checkbox', 'String', 'Remove mean',...
    'Value', removeMean0, 'Units', 'characters', 'Position', [2, 12.5, 40, 1.5], 'HorizontalAlignment', 'left');

% type of detection
detecPan = uipanel(fig, 'Units', 'characters', 'Position', [1, 5.5, 50, 6.5]);
uicontrol(detecPan, 'Style', 'text', 'String', 'Random decrement detection method:',...
    'Units', 'characters', 'Position', [1, 4, 50, 1.7], 'HorizontalAlignment', 'left');
detecMethodNames = {'threshold', 'zero level positive slope', 'zero level negative slope'};
detecMethodInput = uicontrol(detecPan, 'Style', 'popupmenu', 'String', detecMethodNames,...
    'Units', 'characters', 'Position', [2, 2.7, 28, 1.7], 'HorizontalAlignment', 'left');

% threshold
thresholdTxt = uicontrol(detecPan, 'Style', 'text', 'String', 'Threshold:',...
    'Units', 'characters', 'Position', [1, 0.2, 12, 1.4], 'HorizontalAlignment', 'left');
thresholdInput = uicontrol(detecPan, 'Style', 'edit', 'String', thresholdValue0,...
    'Value', true, 'Units', 'characters', 'Position', [13, 0.2, 7, 1.6], 'HorizontalAlignment', 'left');
thresholdMethodNames = {'standard dev', 'absolute'};
thresholdMethodInput = uicontrol(detecPan, 'Style', 'popupmenu', 'String', thresholdMethodNames,...
    'Value', true, 'Units', 'characters', 'Position', [25, 0.3, 18, 1.7], 'HorizontalAlignment', 'left');

% lag time
uicontrol(fig, 'Style', 'text', 'String', 'Lag time:',...
    'Value', true, 'Units', 'characters', 'Position', [2, 3.5, 10, 1.4], 'HorizontalAlignment', 'left');
lagTimeInput = uicontrol(fig, 'Style', 'edit', 'String', lagTime0,...
    'Value', true, 'Units', 'characters', 'Position', [12, 3.5, 8, 1.6], 'HorizontalAlignment', 'left');

% ok button
okInput = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'OK',...
    'Units', 'characters', 'Position', [21, 0.5, 10, 2.2]);

% focus
uicontrol(lagTimeInput);



%% callback fct

KTrdt = [];
lagTime = [];
removeMean = [];
computingRDT = false;

menuCallback;


    function menuCallback(~,~)
        computingRDT = true;
        set(okInput, 'Enable', 'off');
        set(fig, 'pointer', 'watch');
        drawnow
        
        %%%% input
        
        refCh = get(refChInput, 'Value');
        removeMean = get(removeMeanInput, 'Value');
        detecMethod = detecMethodNames{get(detecMethodInput, 'Value')};
        switch detecMethod
            case 'threshold'
                set(thresholdTxt, 'Enable', 'on');
                set(thresholdInput, 'Enable', 'on');
                set(thresholdMethodInput, 'Enable', 'on');
                thresholdValue = str2double(get(thresholdInput, 'String'));
            otherwise
                set(thresholdTxt, 'Enable', 'off');
                set(thresholdInput, 'Enable', 'off');
                set(thresholdMethodInput, 'Enable', 'off');
                thresholdValue = 0;
        end
        thresholdMethod = thresholdMethodNames{get(thresholdMethodInput, 'Value')};
        switch thresholdMethod
            case 'absolute'
                thresholdString = num2str(thresholdValue);
            case 'standard dev'
                switch detecMethod
                    case 'threshold'
                        thresholdString = [num2str(thresholdValue), ' * std dev'];
                    otherwise
                        thresholdString = num2str(thresholdValue);
                end
                thresholdValue = thresholdValue * StdDev(refCh);
        end
        lagTime = str2double(get(lagTimeInput, 'String'));
        
        
        %%%% RDT times
        
        KTrdt = getRDTtimes(t, X - removeMean * mean(X, 2), lagTime, refCh, detecMethod, thresholdValue);
        
        %%%% display update
        
        plotRDTtimes.update(removeMean, refCh, thresholdValue, thresholdString, KTrdt)
        
        %%%% end
        set(fig, 'pointer', 'arrow');
        set(okInput, 'Enable', 'on');
        computingRDT = false;
    end


refChInput.Callback = @menuCallback;
removeMeanInput.Callback = @menuCallback;
detecMethodInput.Callback = @menuCallback;
thresholdInput.Callback = @menuCallback;
thresholdMethodInput.Callback = @menuCallback;
lagTimeInput.Callback = @menuCallback;


%% ok callback

okInput.Callback = @okCallback;

trdt = [];
Xrdt = [];
axRdt = gobjects(1, 0);

    function okCallback(~,~)
        if computingRDT
            return
        end
        
        if lagTime < mean(diff(t))
            errordlg('lag time < dt', 'Error', 'modal');
            uicontrol(lagTimeInput);
            return
        end
        if isempty(KTrdt)
            errordlg('threshold not reached', 'Error', 'modal');
            uicontrol(thresholdInput);
            return
        end
        
        set(okInput, 'Enable', 'off');
        set(fig, 'pointer', 'watch');
        drawnow
        
        % rdt
        [trdt, Xrdt] = getRDT(t, X - removeMean * mean(X, 2), lagTime, KTrdt);
        
        % plot
        axRdt = axes(figure('Name', 'Random decrement signature'));
        plot(axRdt, trdt, Xrdt);
%         set(axRdt, 'XLim', [trdt(1), trdt(end)]);
        legend(axRdt, refChNames);
        xlabel(axRdt, 'Time [s]');
        ylabel(axRdt, ['Random decrement signature [', signalUnit, ']']);
        
        % end
        try
            delete(fig);
        catch
        end
    end



waitfor(fig);



end

