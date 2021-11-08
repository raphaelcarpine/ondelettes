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

Nbuttons = 7;
Hbuttons = 2.5; % characters

Hpans = [Nbuttons*Hbuttons, 4, 4];

fig = figure('Name', 'Linear Filter Menu', 'numbertitle', 'off', 'Resize', 'off');
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

paramPan.Units = 'characters';
Hbuttons = paramPan.Position(4)/Nbuttons;

%% param panel

removeMeanCheck = uicontrol('Parent', paramPan, 'Style', 'checkbox', 'String', 'remove mean',...
    'Units', 'normalized', 'Position', [0.03, 0, 0.8, 0]);
removeMeanCheck.Units = 'characters'; removeMeanCheck.Position([2 4]) = Hbuttons*[6 1];
smoothSignalCheck = uicontrol('Parent', paramPan, 'Style', 'checkbox', 'String', 'smooth signal',...
    'Units', 'normalized', 'Position', [0.03, 0, 0.8, 0]);
smoothSignalCheck.Units = 'characters'; smoothSignalCheck.Position([2 4]) = Hbuttons*[5 1];
removeLocalMeanCheck = uicontrol('Parent', paramPan, 'Style', 'checkbox', 'String', 'remove local mean',...
    'Units', 'normalized', 'Position', [0.03, 0, 0.8, 0]);
removeLocalMeanCheck.Units = 'characters'; removeLocalMeanCheck.Position([2 4]) = Hbuttons*[4 1];
highPassCheck = uicontrol('Parent', paramPan, 'Style', 'checkbox', 'String', 'high pass',...
    'Units', 'normalized', 'Position', [0.03, 0, 0.25, 0]);
highPassCheck.Units = 'characters'; highPassCheck.Position([2 4]) = Hbuttons*[3 1];
lowPassCheck = uicontrol('Parent', paramPan, 'Style', 'checkbox', 'String', 'low pass',...
    'Units', 'normalized', 'Position', [0.03, 0, 0.25, 0]);
lowPassCheck.Units = 'characters'; lowPassCheck.Position([2 4]) = Hbuttons*[2 1];
deriveCheck = uicontrol('Parent', paramPan, 'Style', 'checkbox', 'String', 'derive',...
    'Units', 'normalized', 'Position', [0.03, 0, 0.8, 0]);
deriveCheck.Units = 'characters'; deriveCheck.Position([2 4]) = Hbuttons*[1 1];
integrateCheck = uicontrol('Parent', paramPan, 'Style', 'checkbox', 'String', 'integrate',...
    'Units', 'normalized', 'Position', [0.03, 0, 0.8, 0]);
integrateCheck.Units = 'characters'; integrateCheck.Position([2 4]) = Hbuttons*[0 1];

localMeanWindows = {'rectangular', 'gaussian'};
smoothSignalType = uicontrol('Parent', paramPan, 'Style', 'popupmenu', 'String', localMeanWindows,...
    'Units', 'normalized', 'Position', [0.43, 3/Nbuttons, 0.25, 0.25]);
smoothSignalType.Units = 'characters'; smoothSignalType.Position([2 4]) = Hbuttons*[5 0.8];
localMeanType = uicontrol('Parent', paramPan, 'Style', 'popupmenu', 'String', localMeanWindows,...
    'Units', 'normalized', 'Position', [0.43, 3/Nbuttons, 0.25, 0.25]);
localMeanType.Units = 'characters'; localMeanType.Position([2 4]) = Hbuttons*[4 0.8];
filtersNames = {'butterworth'};
highPassType = uicontrol('Parent', paramPan, 'Style', 'popupmenu', 'String', filtersNames,...
    'Units', 'normalized', 'Position', [0.26, 3/Nbuttons, 0.25, 0.25]);
highPassType.Units = 'characters'; highPassType.Position([2 4]) = Hbuttons*[3 0.8];
lowPassType = uicontrol('Parent', paramPan, 'Style', 'popupmenu', 'String', filtersNames,...
    'Units', 'normalized', 'Position', [0.26, 2/Nbuttons, 0.25, 0.25]);
lowPassType.Units = 'characters'; lowPassType.Position([2 4]) = Hbuttons*[2 0.8];

smoothSignalTimeStr = uicontrol('Parent', paramPan, 'Style', 'text', 'String', 'time:',...
    'Units', 'normalized', 'Position', [0.72, 0.38, 0.1, 0.2]);
smoothSignalTimeStr.Units = 'characters'; smoothSignalTimeStr.Position([2 4]) = Hbuttons*[5 0.7];
smoothSignalTime = uicontrol('Parent', paramPan, 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.82, 0.38, 0.13, 0.23]);
smoothSignalTime.Units = 'characters'; smoothSignalTime.Position([2 4]) = Hbuttons*[5 0.8];
localMeanTimeStr = uicontrol('Parent', paramPan, 'Style', 'text', 'String', 'time:',...
    'Units', 'normalized', 'Position', [0.72, 0.38, 0.1, 0.2]);
localMeanTimeStr.Units = 'characters'; localMeanTimeStr.Position([2 4]) = Hbuttons*[4 0.7];
localMeanTime = uicontrol('Parent', paramPan, 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.82, 0.38, 0.13, 0.23]);
localMeanTime.Units = 'characters'; localMeanTime.Position([2 4]) = Hbuttons*[4 0.8];
highPassFreqStr = uicontrol('Parent', paramPan, 'Style', 'text', 'String', 'freq:',...
    'Units', 'normalized', 'Position', [0.54, 0.38, 0.1, 0.2]);
highPassFreqStr.Units = 'characters'; highPassFreqStr.Position([2 4]) = Hbuttons*[3 0.7];
highPassFreq = uicontrol('Parent', paramPan, 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.64, 0.38, 0.13, 0.23]);
highPassFreq.Units = 'characters'; highPassFreq.Position([2 4]) = Hbuttons*[3 0.8];
lowPassFreqStr = uicontrol('Parent', paramPan, 'Style', 'text', 'String', 'freq:',...
    'Units', 'normalized', 'Position', [0.54, 0.05, 0.1, 0.2]);
lowPassFreqStr.Units = 'characters'; lowPassFreqStr.Position([2 4]) = Hbuttons*[2 0.7];
lowPassFreq = uicontrol('Parent', paramPan, 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.64, 0.05, 0.13, 0.23]);
lowPassFreq.Units = 'characters'; lowPassFreq.Position([2 4]) = Hbuttons*[2 0.8];

highPassNStr = uicontrol('Parent', paramPan, 'Style', 'text', 'String', 'n:',...
    'Units', 'normalized', 'Position', [0.8, 0.38, 0.05, 0.2]);
highPassNStr.Units = 'characters'; highPassNStr.Position([2 4]) = Hbuttons*[3 0.7];
highPassN = uicontrol('Parent', paramPan, 'Style', 'edit', 'String', defaultN,...
    'Units', 'normalized', 'Position', [0.85, 0.38, 0.1, 0.23]);
highPassN.Units = 'characters'; highPassN.Position([2 4]) = Hbuttons*[3 0.8];
lowPassNStr = uicontrol('Parent', paramPan, 'Style', 'text', 'String', 'n:',...
    'Units', 'normalized', 'Position', [0.8, 0.05, 0.05, 0.2]);
lowPassNStr.Units = 'characters'; lowPassNStr.Position([2 4]) = Hbuttons*[2 0.7];
lowPassN = uicontrol('Parent', paramPan, 'Style', 'edit', 'String', defaultN,...
    'Units', 'normalized', 'Position', [0.85, 0.05, 0.1, 0.23]);
lowPassN.Units = 'characters'; lowPassN.Position([2 4]) = Hbuttons*[2 0.8];


    function removeMean()
        for k = 1:N
            X{k} = X{k} - mean(X{k});
        end
    end


    function removeLocalMean()
        timeMean = eval(get(localMeanTime, 'String'));
        windowType = localMeanWindows{get(localMeanType, 'Value')};
        for k = 1:N
            X{k} = X{k} - getSmoothSignal(T0{k}, X{k}, windowType, timeMean);
        end
    end


    function smoothSignal()
        timeMean = eval(get(smoothSignalTime, 'String'));
        windowType = localMeanWindows{get(smoothSignalType, 'Value')};
        for k = 1:N
            X{k} = getSmoothSignal(T0{k}, X{k}, windowType, timeMean);
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

    function applyDerive()
        for k = 1:N
            dt = mean(diff(T0{k}));
            if max(abs(diff(T0{k})/dt -1)) > 1e-4
                error('sampling problem');
            end
            
            dX = 1/dt * diff(X{k});
            X{k} = [dX(1), 1/2*(dX(1:end-1) + dX(2:end)), dX(end)];
        end
    end

    function applyIntegrate()
        for k = 1:N
            dt = mean(diff(T0{k}));
            if max(abs(diff(T0{k})/dt -1)) > 1e-4
                error('sampling problem');
            end
            
            iX = zeros(size(X{k}));
            for i = 2:length(X{k})
                iX(i) = iX(i-1) + dt * X{k}(i-1);
            end
            X{k} = iX;
        end
    end



%% filter panel

buttonFilter = uicontrol('Parent',filterPan, 'Units', 'normalized','Style','pushbutton',...
    'String', 'filter', 'Position', [0.25, 0.06, 0.5, 0.88]);

    function applyFilters()
        if get(removeMeanCheck, 'Value')
            removeMean();
        end
        if get(smoothSignalCheck, 'Value')
            smoothSignal();
        end
        if get(removeLocalMeanCheck, 'Value')
            removeLocalMean();
        end
        if get(highPassCheck, 'Value')
            applyHighPassFilter();
        end
        if get(lowPassCheck, 'Value')
            applyLowPassFilter();
        end
        if get(deriveCheck, 'Value')
            applyDerive();
        end
        if get(integrateCheck, 'Value')
            applyIntegrate();
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

