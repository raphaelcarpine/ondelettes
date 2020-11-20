function ax = plotMeanAndShocks(t, meanWvlt, meanWvlttEdgeEffects, shockIndexes, threshold, meanT, meanScale, figName)
%PLOTMEANANDSHOCKS Summary of this function goes here
%   Detailed explanation goes here

%% mean plot

fig = figure('Name', figName);
ax = axes(fig);
hold(ax, 'on');
plt0 = plot(ax, t, meanWvlttEdgeEffects, '--');
plt = plot(ax, t, meanWvlt, 'Color', get(plt0, 'Color'));
set(ax, 'YScale', meanScale);
xlabel(ax, 'Time [s]');
ylabel(ax, 'Average F(CWT(t,~))');


%% shocks

shockTimes = t(shockIndexes);
shockValues = meanWvlt(shockIndexes);
shockTimeNames = ['t_{'.*ones(length(shockTimes), 1), num2str(transpose(1:length(shockTimes))), '}'.*ones(length(shockTimes), 1)];
shockTimeNames = cellstr(shockTimeNames);

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
plot(ax, [t(1), t(end)], [threshold, threshold], '-', 'Color', colorThreshold);
text(ax, t(1), threshold, '  threshold', 'Color', colorThreshold, 'VerticalAlignment', 'bottom');

% mean
colorMean = 0*[1 1 1];
plot(ax, [t(1), t(end)], [meanT, meanT], '--', 'Color', colorMean);
text(ax, t(1), meanT, '  mean', 'Color', colorMean, 'VerticalAlignment','bottom');

end

