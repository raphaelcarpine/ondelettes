function [Xav, Yav, stdY, K, Xlims] = ScatterAveragingMenu(varargin)
%WaveletMenu Summary of this function goes here
%   Detailed explanation goes here

figureName = 'Scatter Averaging';

p = inputParser;

defaultFunction = '@(x, y) averagingScatter2(x, y, 1000)';
defaultXaxis = '@(x) x';
defaultYaxis = '@(y) y';

addOptional(p, 'Function', defaultFunction);
addOptional(p, 'Xaxis', defaultXaxis);
addOptional(p, 'Yaxis', defaultYaxis);

parse(p, varargin{:});

func = p.Results.Function;
Xaxis = p.Results.Xaxis;
Yaxis = p.Results.Yaxis;


%% fig

% lignes pans
Hpans = [4, 7, 8.5, 3];

fig = figure('Name', figureName, 'numbertitle', 'off', 'Resize', 'off');
fig.Units = 'characters';
fig.Position(3) = 70;
fig.Position(4) = sum(Hpans);
fig.Position(2) = fig.Position(2) - fig.Position(4) - 2;
fig.MenuBar = 'none';



%% bouton regression et panneaux param et sorties

eqPan = uipanel('Parent',fig, 'Units', 'normalized');
linePan = uipanel('Parent',fig, 'Units', 'normalized');
optionsPan = uipanel('Parent',fig, 'Units', 'normalized', 'Title', 'Plot');
buttonAverage = uicontrol('Parent',fig, 'Units', 'normalized','Style','pushbutton',...
    'String', 'Average');

margin = 0.02;
hpans = (1 - margin*(length(Hpans)+1)) * Hpans/sum(Hpans);
zpans = margin * ones(size(hpans));
for iz = length(zpans)-1:-1:1
    zpans(iz) = zpans(iz+1) + hpans(iz+1) + margin;
end

linePan.Position = [0.02, zpans(1), 0.96, hpans(1)];
eqPan.Position = [0.02, zpans(2), 0.96, hpans(2)];
optionsPan.Position = [0.02, zpans(3), 0.96, hpans(3)];
buttonAverage.Position = [0.02, zpans(4), 0.96, hpans(4)];



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
prevBut.Position = [0.6, 0.1, 0.15, 0.8];
nextBut.Position = [0.75, 0.1, 0.15, 0.8];


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
            lines = findobj({'Type', 'line'},'-or', {'Type', 'scatter'});
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



%% function

funcStr = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','text', 'String', 'averaging function:');
funcEdit = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','edit', 'String', func);
XaxisStr = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','text', 'String', 'x-axis:');
XaxisEdit = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','edit', 'String', Xaxis);
YaxisStr = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','text', 'String', 'y-axis:');
YaxisEdit = uicontrol('Parent', eqPan, 'Units', 'normalized','Style','edit', 'String', Yaxis);

nLign = 3;
marge = 0.02;
h = (1-(nLign+1)*marge)/nLign;
H = h+marge;
funcStr.Position = [0.01, 1-H, 0.2, h];
funcEdit.Position = [0.22, 1-H, 0.77, h];
XaxisStr.Position = [0.01, 1-2*H, 0.2, h];
XaxisEdit.Position = [0.22, 1-2*H, 0.77, h];
YaxisStr.Position = [0.01, 1-3*H, 0.2, h];
YaxisEdit.Position = [0.22, 1-3*H, 0.77, h];


%% options

options = {'onaxes', 'std', 'xlims', 'K'};
nopt = length(options);

optionsStr = struct;
optionsStr.onaxes = 'current axes';
optionsStr.std = 'standard deviation';
optionsStr.xlims = 'samples limits';
optionsStr.K = 'samples sizes';

optBut = struct;
for kopt = 1:length(options)
    opt = options{kopt};
    optBut.(opt) = uicontrol('Parent', optionsPan, 'Units', 'normalized','Style','checkbox',...
        'String', optionsStr.(opt), 'Position', [0.01, (nopt-kopt)/nopt, 0.98, 1/nopt]);
end
optBut.onaxes.Value = true;
optBut.std.Value = true;
optBut.xlims.Value = false;
optBut.K.Value = false;

% delete plots button
onAxesPlots = [];
deleteBut = uicontrol('Parent', optionsPan, 'Units', 'normalized','Style','pushbutton',...
    'String', 'delete plots', 'Position', [0.7, (nopt-1)/nopt, 0.29, 1/nopt]);

    function deletePlots()
        delete(onAxesPlots);
        onAxesPlots = [];
    end

deleteBut.Callback = @(~,~) deletePlots();




%% donnees

X = nan;
Y = nan;
ax = 0;

    function ok = updateXYAxes()
        %         line = findobj('Type', 'line');
        if ~isempty(line)
            %             line = line(1);
            X = get(line, 'XData');
            Y = get(line, 'YData');
            
            X = X(~isnan(Y)); % on enlève les valeurs inappropriées
            Y = Y(~isnan(Y));
            X = X(~isnan(X));
            Y = Y(~isnan(X));
            
            ax = get(line, 'Parent');
            ok = ~isempty(X);
        else
            ok = false;
        end
    end



%% averaging

Xav = nan;
Yav = nan;
stdY = nan;
K = nan;
Xlims = nan;

    function computeAverage()
        if ~updateXYAxes()
            warning('no line selected');
            return;
        end
        
        % upload
        func = str2func(funcEdit.String);
        Xaxis = str2func(XaxisEdit.String);
        Yaxis = str2func(YaxisEdit.String);
        plotSTD = optBut.std.Value;
        plotXlims = optBut.xlims.Value;
        
        %
        set(YaxisEdit, 'ForegroundColor', [0.5 0.5 0.5]);
        set(fig, 'pointer', 'watch');
        drawnow;
        
        % averaging
        [Xav, Yav, stdY, K, Xlims] = func(Xaxis(X), Yaxis(Y));
        
        % end
        set(YaxisEdit, 'ForegroundColor', [0 0 0]);
        set(fig, 'pointer', 'arrow');
        drawnow
        
        plotAverage(Xav, Yav, stdY, Xlims, plotSTD, plotXlims);
        
        if optBut.K.Value
            plotSampleSize(Xav, K);
        end
        
        lineSelect.Value = false;
        selectFunction(lineSelect.Value);
    end


    function plotAverage(Xav, Yav, stdY, Xlims, plotSTD, plotXlims)
        
        if optBut.onaxes.Value
            plotAxes = ax;
            hold(plotAxes, 'on');
            ylim(ax, 'manual');
            
            onAxesPlots = [onAxesPlots, plot(plotAxes, Xav, Yav, 'r+', 'LineWidth', 2)];
            
            if plotSTD
                for indXav = 1:length(Xav)
                    onAxesPlots = [onAxesPlots,...
                        plot(plotAxes, [Xav(indXav), Xav(indXav)], Yav(indXav) + stdY(indXav) * [-1 1], '-r_')];
                end
            end
            
            if plotXlims
                for indXlims = 1:length(Xlims)
                    onAxesPlots = [onAxesPlots, xline(plotAxes, Xlims(indXlims), '--')];
                end
            end
            
            %uistack(onAxesPlots(end), 'bottom'); %ligne derrière/devant
            hold(plotAxes, 'off');
        else
            plotAxes = axes(figure);
            hold(plotAxes, 'on');
            
            plot(plotAxes, X, Y, '*');
            plot(plotAxes, Xav, Yav, 'r+', 'LineWidth', 2);
            
            if plotSTD
                for indXav = 1:length(Xav)
                    plot(plotAxes, [Xav(indXav), Xav(indXav)], Yav(indXav) + stdY(indXav) * [-1 1], '-r_');
                end
            end
            
            if plotXlims
                for indXlims = 1:length(Xlims)
                    xline(plotAxes, Xlims(indXlims), '--');
                end
            end
            
            hold(plotAxes, 'off');
        end
    end



buttonAverage.Callback = @(~,~) computeAverage();

%% sample size

    function plotSampleSize(Xav, K)
        axBar = axes(figure);
        bar(axBar, Xav, K);
        xlabel(axBar, 'x');
        ylabel(axBar, 'sample size');
    end


%% raccourcis

    function keyboardShortcuts(~, event)
        if isequal(event.Key, 'return') || isequal(event.Key, 'space')
            computeAverage('reg');
        elseif isequal(event.Key, 'shift') || isequal(event.Key, 's')
            lineSelect.Value = ~ lineSelect.Value;
            selectFunction(lineSelect.Value);
        elseif isequal(event.Key, 'p')
            plotFunc();
        elseif isequal(event.Key, 'd')
            deletePlots();
        elseif isequal(event.Key, 'leftarrow')
            nextprev(-1);
        elseif isequal(event.Key, 'rightarrow')
            nextprev(1);
        elseif isequal(event.Key, 'delete') || isequal(event.Key, 'backspace')
            close(fig);
        end
    end
set(fig, 'KeyPressFcn', @keyboardShortcuts);


%% fermeture

if nargout > 0
    waitfor(fig);
end


end