function ax = plotMeanAndShocks(t, meanWvlt, meanWvlttEdgeEffects, shockIndexes, threshold, meanT, meanScale, deltaT, figName)
%PLOTMEANANDSHOCKS Summary of this function goes here
%   Detailed explanation goes here

displayDeltaTfmin = true;
displayDeltaTfmax = true;

%% mean plot

fig = figure('Name', figName);
ax = axes(fig);
hold(ax, 'on');
plt0 = plot(ax, t, meanWvlttEdgeEffects, '--');
plt = plot(ax, t, meanWvlt, 'Color', get(plt0, 'Color'));
set(ax, 'YScale', meanScale);
set(ax, 'XLim', [t(1), t(end)]);
xlabel(ax, 'Time [s]');
ylabel(ax, 'Average F(CWT(t,~))');


%% shocks

shockTimes = t(shockIndexes);
shockValues = meanWvlt(shockIndexes);
shockTimeNames = cell(length(shockTimes), 1);
for k_shock = 1:length(shockTimes)
    shockTimeNames{k_shock} = ['t_{', num2str(k_shock), '}'];
end

% graph
plot(ax, shockTimes, meanWvlt(shockIndexes), '+',...
    'Color', get(plt, 'Color'), 'LineWidth', 2*get(plt, 'LineWidth'));
for kt = 1:length(shockTimes)
    plot(ax, [shockTimes(kt), shockTimes(kt)], [0, shockValues(kt)], ':', 'Color', 0.5*[1 1 1]);
end

% ticks
% XTick2 = get(ax, 'XTick');
% XTicklabel2 = get(ax, 'XTickLabel');
% XTicklabel2 = XTicklabel2(~any(XTick2 == shockTimes'));
% XTick2 = XTick2(~any(XTick2 == shockTimes'));
% XTick2 = [XTick2, shockTimes];
% XTicklabel2 = [XTicklabel2; shockTimeNames];
% [XTick2, Ixtick] = sort(XTick2);
% XTicklabel2 = XTicklabel2(Ixtick);
% set(ax, 'XTick', XTick2);
% set(ax, 'XTickLabel', XTicklabel2);

% times
for k_t = 1:length(shockTimes)
    text(ax, shockTimes(k_t), shockValues(k_t), shockTimeNames{k_t}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end


%% mean & threshold

% threshold
colorThreshold = 'r';
yline(ax, threshold, '-', 'threshold', 'Color', colorThreshold, 'LabelHorizontalAlignment', 'left');

% mean
colorMean = 0*[1 1 1];
yline(ax, meanT, '--', 'average', 'Color', colorMean, 'LabelHorizontalAlignment', 'right');


%% delta T
deltaTfminTxt = text(ax, nan, nan, [' \Deltat', sprintf(' for f_{min} = %.2gs ', deltaT(1))],...
    'EdgeColor', [0 0 0], 'BackgroundColor', [1 1 1], 'Margin', 1, 'Visible', displayDeltaTfmin);
deltaTfmaxTxt = text(ax, nan, nan, [' \Deltat', sprintf(' for f_{max} = %.2gs ', deltaT(2))],...
    'EdgeColor', [0 0 0], 'BackgroundColor', [1 1 1], 'Margin', 1, 'Visible', displayDeltaTfmax);
deltaTfminPlt = plot(ax, nan, nan, '-+', 'Color', [0 0 0], 'Visible', displayDeltaTfmin);%,...
%     'LineWidth', 2);
deltaTfmaxPlt = plot(ax, nan, nan, '-+', 'Color', [0 0 0], 'Visible', displayDeltaTfmax);%,...
%     'LineWidth', 2);



    function updateDeltaT(~,~)
        set(ax, 'Units', 'characters');
        sizeChar = ax.InnerPosition(3:4);
        set(ax, 'Units', 'normalized');
        
        set(deltaTfminTxt, 'Units', 'characters');
        set(deltaTfmaxTxt, 'Units', 'characters');
        deltaTfminTxt.Position(1) = 3;
        deltaTfminTxt.Position(2) = sizeChar(2) - 3;
        deltaTfmaxTxt.Position(1) = 3;
        deltaTfmaxTxt.Position(2) = sizeChar(2) - 6;
        set(deltaTfminTxt, 'Units', 'data');
        set(deltaTfmaxTxt, 'Units', 'data');
        
        set(deltaTfminPlt, 'XData', deltaTfminTxt.Position(1) + [0, deltaT(1)],...
            'YData', deltaTfminTxt.Position(2) * [1 1]);
        set(deltaTfmaxPlt, 'XData', deltaTfmaxTxt.Position(1) + [0, deltaT(2)],...
            'YData', deltaTfmaxTxt.Position(2) * [1 1]);
        
        set(deltaTfminTxt, 'Units', 'characters');
        set(deltaTfmaxTxt, 'Units', 'characters');
        deltaTfminTxt.Position(1) = 2;
        deltaTfminTxt.Position(2) = sizeChar(2) - 1.5;
        deltaTfmaxTxt.Position(1) = 2;
        deltaTfmaxTxt.Position(2) = sizeChar(2) - 4.5;
    end

updateDeltaT;

addlistener(ax, 'XLim', 'PostSet', @updateDeltaT);
addlistener(ax, 'YLim', 'PostSet', @updateDeltaT);
addlistener(ax, 'YScale', 'PostSet', @updateDeltaT);
fig.SizeChangedFcn = @updateDeltaT;

%% options

optionsMenu = uimenu(fig, 'Text', 'OPTIONS');

%%%% sclaes
% amplitudescale
amplMenu = uimenu(optionsMenu, 'Text', 'Amplitude scale');
linAmplMenu = uimenu(amplMenu, 'Text', 'lin');
logAmplMenu = uimenu(amplMenu, 'Text', 'log');
switch meanScale
    case 'lin'
        set(linAmplMenu, 'Checked', 'on');
    case 'log'
        set(logAmplMenu, 'Checked', 'on');
end

    function linAmplMenuCallback(~, ~)
        set(linAmplMenu, 'Checked', 'on');
        set(logAmplMenu, 'Checked', 'off');
        set(ax, 'YScale', 'lin');
    end
    function logAmplMenuCallback(~, ~)
        set(linAmplMenu, 'Checked', 'off');
        set(logAmplMenu, 'Checked', 'on');
        set(ax, 'YScale', 'log');
    end
linAmplMenu.Callback = @linAmplMenuCallback;
logAmplMenu.Callback = @logAmplMenuCallback;


%%%% Delta T display
deltaTfminMenu = uimenu(optionsMenu, 'Text', 'Display Delta t for fmin', 'Checked', displayDeltaTfmin,...
    'Separator', 'on');
deltaTfmaxMenu = uimenu(optionsMenu, 'Text', 'Display Delta t for fmax', 'Checked', displayDeltaTfmax);

    function deltaTfminCallback(~, ~)
        if displayDeltaTfmin
            set(deltaTfminMenu, 'Checked', 'off');
        else
            set(deltaTfminMenu, 'Checked', 'on');
        end
        displayDeltaTfmin = ~displayDeltaTfmin;
        set(deltaTfminPlt, 'Visible', displayDeltaTfmin);
        set(deltaTfminTxt, 'Visible', displayDeltaTfmin);
    end

    function deltaTfmaxCallback(~, ~)
        if displayDeltaTfmax
            set(deltaTfmaxMenu, 'Checked', 'off');
        else
            set(deltaTfmaxMenu, 'Checked', 'on');
        end
        displayDeltaTfmax = ~displayDeltaTfmax;
        set(deltaTfmaxPlt, 'Visible', displayDeltaTfmax);
        set(deltaTfmaxTxt, 'Visible', displayDeltaTfmax);
    end

deltaTfminMenu.Callback = @deltaTfminCallback;
deltaTfmaxMenu.Callback = @deltaTfmaxCallback;

end

