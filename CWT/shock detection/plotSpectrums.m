function plotSpectrums(Freqs, Spectrums, spectrumFrequencyScale, spectrumScale, figName,...
    plotAverageShockSpectrum, plotAverageSpectrum, plotUnderThresholdAverage, plotAboveThresholdAverage, customLegend)
%PLOTSPECTRUMS Summary of this function goes here
%   Detailed explanation goes here

if nargin <=5
    plotAverageShockSpectrum = false;
    plotAverageSpectrum = false;
    plotUnderThresholdAverage = false;
    plotAboveThresholdAverage = false;
end
if nargin <= 9
    customLegend = {};
end

fig = figure('Name', figName);
ax = axes(fig);
hold(ax, 'on');
set(ax, 'XScale', spectrumFrequencyScale, 'YScale', spectrumScale);
set(ax, 'XLim', [Freqs(1), Freqs(end)]);
xlabel(ax, 'Frequency [Hz]');
ylabel(ax, 'Squared amplitude CWT');


%% plot

if plotAverageShockSpectrum || plotAverageSpectrum || plotUnderThresholdAverage || plotAboveThresholdAverage
    plotIndexes = [plotAverageShockSpectrum, plotAverageSpectrum, plotUnderThresholdAverage, plotAboveThresholdAverage];
    plotIndexes = cumsum(plotIndexes) .* plotIndexes;
    
    % legend
    legendStrings = {'average shocks spectrum', 'average spectrum', 'average under threshold spectrum', 'average above threshold spectrum'};
    legendStrings = legendStrings(logical(plotIndexes));
    
    % plot
    plts = gobjects(0, 1);
    for k = 1:length(plotIndexes)
        if plotIndexes(k) == 0
            continue
        end
        plts = [plts; plot(ax, Freqs, Spectrums(plotIndexes(k), :), 'LineWidth', 2)];
        
        if k == 2
            set(plts(end), 'LineStyle','--');
        elseif k == 3 || k == 4
            set(plts(end), 'LineStyle',':');
        end
        
    end
    
elseif isempty(customLegend)
    % legend
    shockTimeNames = cell(size(Spectrums, 1), 1);
    for k_shock = 1:size(Spectrums, 1)
        shockTimeNames{k_shock} = ['t_{', num2str(k_shock), '}'];
    end
    legendStrings = shockTimeNames;
    
    plts = plot(ax, Freqs, Spectrums);
    
else
    legendStrings = customLegend;
    plts = plot(ax, Freqs, Spectrums);
end

% x-axis
if any(Spectrums < 0, 'all')
    yline(ax, 0, '--');
end



% legend
maxLegendLines = 10;

legendOn = length(legendStrings) <= maxLegendLines;
if legendOn
    legend(ax, legendStrings, 'Location', 'best');
end

additionalLegend = gobjects(1, 0);

% data tips
for k_plot = 1:length(plts)
    plt = plts(k_plot);
    plt.DataTipTemplate.DataTipRows = [plt.DataTipTemplate.DataTipRows(1); plt.DataTipTemplate.DataTipRows];
    plt.DataTipTemplate.DataTipRows(1).Label = legendStrings{k_plot};
    plt.DataTipTemplate.DataTipRows(1).Value = '';
    plt.DataTipTemplate.DataTipRows(2).Label = 'freq';
    plt.DataTipTemplate.DataTipRows(3).Label = 'ampl';
end



%% options menu

optionsMenu = uimenu(fig, 'Text', 'OPTIONS');

%%%% sclaes
% amplitudescale
amplMenu = uimenu(optionsMenu, 'Text', 'Amplitude scale');
linAmplMenu = uimenu(amplMenu, 'Text', 'lin');
logAmplMenu = uimenu(amplMenu, 'Text', 'log');
switch spectrumScale
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

% frequency scale
freqMenu = uimenu(optionsMenu, 'Text', 'Frequency scale');
linFreqMenu = uimenu(freqMenu, 'Text', 'lin');
logFreqMenu = uimenu(freqMenu, 'Text', 'log');
switch spectrumFrequencyScale
    case 'lin'
        set(linFreqMenu, 'Checked', 'on');
    case 'log'
        set(logFreqMenu, 'Checked', 'on');
end

    function linFreqMenuCallback(~, ~)
        set(linFreqMenu, 'Checked', 'on');
        set(logFreqMenu, 'Checked', 'off');
        set(ax, 'XScale', 'lin');
    end
    function logFreqMenuCallback(~, ~)
        set(linFreqMenu, 'Checked', 'off');
        set(logFreqMenu, 'Checked', 'on');
        set(ax, 'XScale', 'log');
    end
linFreqMenu.Callback = @linFreqMenuCallback;
logFreqMenu.Callback = @logFreqMenuCallback;


%%%% spectrum selection
showAllMenu = uimenu(optionsMenu, 'Text', 'Show all spectrums', 'Checked', 'on',  'Separator', 'on');

selectSpectrumMenu = uimenu(optionsMenu, 'Text', 'Selected displayed spectrums');

    function showAllMenuCallback(~, ~)
        set(showAllMenu, 'Checked', 'on');
        set(selectSpectrumMenu, 'Checked', 'off');
        set(plts, 'Visible', 'on')
        delete(additionalLegend);
    end

    function selectSpectrumMenuCallback(~, ~)
        % dialog box
        validAnswer = false;
        while ~validAnswer
            dlgtitle = 'Spectrum selection';
            prompt = {'Only display spectrums #:'};
            dims = [1 35];
            definput = {''};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            if isempty(answer)
                return
            end
            
            visibleSpectrums = str2num(answer{1});
            if isempty(visibleSpectrums) || ~all(mod(visibleSpectrums, 1) == 0)...
                    || min(visibleSpectrums) < 1 || max(visibleSpectrums) > length(plts)
                errorFig = errordlg('Incorrect input format', 'Error');
                waitfor(errorFig);
            else
                validAnswer = true;
            end
        end
        
        
        set(showAllMenu, 'Checked', 'off');
        set(selectSpectrumMenu, 'Checked', 'on');
        
        set(plts, 'Visible', 'off');
        set(plts(visibleSpectrums), 'Visible', 'on')
        
        delete(additionalLegend);
        if ~legendOn && sum(logical(visibleSpectrums)) <= maxLegendLines
            additionalLegend = legend(ax, plts(visibleSpectrums), legendStrings(visibleSpectrums), 'Location', 'best');
        end
    end

showAllMenu.Callback = @showAllMenuCallback;
selectSpectrumMenu.Callback = @selectSpectrumMenuCallback;

end

