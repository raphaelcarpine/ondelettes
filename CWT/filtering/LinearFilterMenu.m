function LinearFilterMenu(dataPlot)
%LINEARFILTERMENU Summary of this function goes here
%   dataPlot = plot(t, x)

%%
if nargin == 0
    figure;
    t = 0.01*(0:10000-1);
    dataPlot = plot(t, 0*randn(2, length(t)) + sin(2*pi*[1;10]*t) + [0; 1]);
end

defaultN = 5;

%% initialisation

try
    T0 = get(dataPlot, 'XData');
    X0 = get(dataPlot, 'YData');
    if ~iscell(T0)
        T0 = {T0};
        X0 = {X0};
    end
catch
    warning('no plot selected');
    return
end
N = length(X0);
X = X0;

%% display

Hpans = [10, 4, 4];

fig = figure('Name', 'Linear Filter Menu');
fig.Units = 'characters';
fig.Position(3) = 70;
fig.Position(4) = sum(Hpans);
fig.MenuBar = 'none';


paramPan = uipanel('Parent',fig, 'Units', 'normalized');
filterPan = uipanel('Parent',fig, 'Units', 'normalized');
finalPan = uipanel('Parent',fig, 'Units', 'normalized');

margin = 0.02;
hpans = (1 - margin*(length(Hpans)+1)) * Hpans/sum(Hpans);
zpans = margin * ones(size(hpans));
for iz = length(zpans)-1:-1:1
    zpans(iz) = zpans(iz+1) + hpans(iz+1) + margin;
end

paramPan.Position = [0.02, zpans(1), 0.96, hpans(1)];
filterPan.Position = [0.02, zpans(2), 0.96, hpans(2)];
finalPan.Position = [0.02, zpans(3), 0.96, hpans(3)];

%% param panel

removeMeanCheck = uicontrol('Parent', paramPan, 'Style', 'checkbox', 'String', 'remove mean',...
    'Units', 'normalized', 'Position', [0.03, 0.66, 0.8, 0.33]);
highPassCheck = uicontrol('Parent', paramPan, 'Style', 'checkbox', 'String', 'high pass',...
    'Units', 'normalized', 'Position', [0.03, 0.33, 0.25, 0.33]);
lowPassCheck = uicontrol('Parent', paramPan, 'Style', 'checkbox', 'String', 'low pass',...
    'Units', 'normalized', 'Position', [0.03, 0., 0.25, 0.33]);

ordersNames = {'butterworth'};
highPassType = uicontrol('Parent', paramPan, 'Style', 'popupmenu', 'String', ordersNames,...
    'Units', 'normalized', 'Position', [0.29, 0.33, 0.21, 0.25]);
lowPassType = uicontrol('Parent', paramPan, 'Style', 'popupmenu', 'String', ordersNames,...
    'Units', 'normalized', 'Position', [0.29, 0., 0.21, 0.25]);

uicontrol('Parent', paramPan, 'Style', 'text', 'String', 'freq:',...
    'Units', 'normalized', 'Position', [0.54, 0.38, 0.1, 0.2]);
highPassFreq = uicontrol('Parent', paramPan, 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.64, 0.38, 0.13, 0.23]);
uicontrol('Parent', paramPan, 'Style', 'text', 'String', 'freq:',...
    'Units', 'normalized', 'Position', [0.54, 0.05, 0.1, 0.2]);
lowPassFreq = uicontrol('Parent', paramPan, 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.64, 0.05, 0.13, 0.23]);

uicontrol('Parent', paramPan, 'Style', 'text', 'String', 'n:',...
    'Units', 'normalized', 'Position', [0.8, 0.38, 0.05, 0.2]);
highPassN = uicontrol('Parent', paramPan, 'Style', 'edit', 'String', defaultN,...
    'Units', 'normalized', 'Position', [0.85, 0.38, 0.1, 0.23]);
uicontrol('Parent', paramPan, 'Style', 'text', 'String', 'n:',...
    'Units', 'normalized', 'Position', [0.8, 0.05, 0.05, 0.2]);
lowPassN = uicontrol('Parent', paramPan, 'Style', 'edit', 'String', defaultN,...
    'Units', 'normalized', 'Position', [0.85, 0.05, 0.1, 0.23]);


    function removeMean()
        for k = 1:N
            X{k} = X{k} - mean(X{k});
        end
    end




    function applyHighPassFilter()
        f = eval(get(highPassFreq, 'String'));
        n = eval(get(highPassN, 'String'));
        if isnan(f) || isnan(n)
            error('invalid syntax');
        end
        
        filterType = get(highPassType, 'Value');
        if filterType == 1 %butterworth
            for k = 1:N
                X{k} = butterworthFilter(T0{k}, X{k}, f, 'high', n);
            end
        end
    end

    function applyLowPassFilter()
        f = eval(get(lowPassFreq, 'String'));
        n = eval(get(lowPassN, 'String'));
        if isnan(f) || isnan(n)
            error('invalid syntax');
        end
        
        filterType = get(lowPassType, 'Value');
        if filterType == 1 %butterworth
            for k = 1:N
                X{k} = butterworthFilter(T0{k}, X{k}, f, 'low', n);
            end
        end
    end




%% filter panel

buttonFilter = uicontrol('Parent',filterPan, 'Units', 'normalized','Style','pushbutton',...
    'String', 'filter', 'Position', [0.25, 0.06, 0.5, 0.88]);

    function applyFilters()
        if get(removeMeanCheck, 'Value')
            removeMean();
        end
        if get(highPassCheck, 'Value')
            applyHighPassFilter();
        end
        if get(lowPassCheck, 'Value')
            applyLowPassFilter();
        end
        
        for k = 1:N
            set(dataPlot(k), 'YData', X{k});
        end
        
        set(buttonOK, 'Enable', 'on');
    end
buttonFilter.Callback = @(~,~) applyFilters();




%% final panel

buttonOK = uicontrol('Parent', finalPan, 'Units', 'normalized','Style','pushbutton',...
    'String', 'ok', 'Position', [0.02, 0.06, 0.47, 0.88], 'Enable', 'off');
buttonCancel = uicontrol('Parent', finalPan, 'Units', 'normalized','Style','pushbutton',...
    'String', 'cancel', 'Position', [0.51, 0.06, 0.47, 0.88]);

    function cancelFiltering()
        for k = 1:N
            set(dataPlot(k), 'YData', X0{k});
        end
        X = X0;
    end
buttonCancel.Callback = @(~,~) cancelFiltering();

buttonOK.Callback = @(~,~) close(fig);


%%



end

