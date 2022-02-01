function getPCA(Spectrums, Npc, logScale, scaleByStd, Freqs, spectrumFrequencyScale, plotScatter, plotPC, plotDistribution, varargin)
%GETPCA Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
addParameter(p, 'plotTitleSuffix', '');
addParameter(p, 'QSpectrum', 0);
parse(p, varargin{:});
plotTitleSuffix = p.Results.plotTitleSuffix;
QSpectrum = p.Results.QSpectrum;
if ~isempty(plotTitleSuffix)
    plotTitleSuffix = [' ; ', plotTitleSuffix];
end


spectrumScale = 'lin';

if nargin == 0 % test
    Spectrums = rand(1000, 1) * ( 2 +sin(linspace(0, 3*pi, 200))); % 1000 shocks, 200 freqs
    Spectrums = Spectrums + (rand(1000, 1) < 0.5) * exp(-linspace(-10, 5, 200).^2);
    Spectrums = Spectrums + 0.05 * randn(size(Spectrums));
    Npc = 2;
    scaleByStd = false;
    Freqs = linspace(1, 10, 200);
    spectrumFrequencyScale = 'lin';
    plotScatter = true;
    plotPC = true;
    plotDistribution = true;
end

%% SVD
X = Spectrums;
if logScale
    X = log(X);
end
X = X - mean(X, 1);
if scaleByStd
    X = X ./ std(X);
end

[U, S, W] = svd(X);
W = W .* (-1 + 2*(max(W, [], 1) >= -min(W, [], 1)));
T = X * W;
L = W'*X'*X*W;
L = diag(L);


%% principal components distribution plot

if plotScatter
    if Npc <= 3
        fig = figure('Name', ['PCA projection', plotTitleSuffix]);
        ax = axes(fig);
        %     hold(ax, 'on');
        
        shockTimeNames = cell(size(Spectrums, 1), 1);
        for k_shock = 1:size(Spectrums, 1)
            shockTimeNames{k_shock} = sprintf(' t_{%u}', k_shock);
        end
        
        if Npc == 1
            shockPoints = scatter(ax, T(:, 1), zeros(size(T(:, 1))), 'x', 'LineWidth', 2);
            shockTexts = text(ax, T(:, 1), zeros(size(T(:, 1))), shockTimeNames, 'VerticalAlignment', 'baseline');
            hold(ax, 'on');
            yline(ax, 0, '--');
            xlabel(ax, 'principal component 1');
            ax.YAxis.Visible = 'off';
        elseif Npc == 2
            shockPoints = scatter(ax, T(:, 1), T(:, 2), 'x', 'LineWidth', 2);
            shockTexts = text(ax, T(:, 1), T(:, 2), shockTimeNames, 'VerticalAlignment', 'baseline');
            hold(ax, 'on');
            xline(ax, 0, '--');
            yline(ax, 0, '--');
            xlabel(ax, 'principal component 1');
            ylabel(ax, 'principal component 2');
        elseif Npc == 3
            shockPoints = scatter3(ax, T(:, 1), T(:, 2), T(:, 3), 'x', 'LineWidth', 2);
            shockTexts = text(ax, T(:, 1), T(:, 2), T(:, 3), shockTimeNames, 'VerticalAlignment', 'baseline');
            hold(ax, 'on');
            xline(ax, 0, '--');
            yline(ax, 0, '--');
            plot3([0 0], [0 0], zlim, '--', 'Color', [0 0 0]);
            xlabel(ax, 'principal component 1');
            ylabel(ax, 'principal component 2');
            zlabel(ax, 'principal component 3');
        end
        
        %%%% options menu
        optionsMenu = uimenu(fig, 'Text', 'OPTIONS');
        
        % spectrum selection
        showAllMenu = uimenu(optionsMenu, 'Text', 'Show all shocks', 'Checked', 'on');
        selectShocksMenu = uimenu(optionsMenu, 'Text', 'Selected displayed shocks');
        
        showAllMenu.Callback = @showAllMenuCallback;
        selectShocksMenu.Callback = @selectShocksMenuCallback;
        
    end
end

    function showAllMenuCallback(~, ~)
        set(showAllMenu, 'Checked', 'on');
        set(selectShocksMenu, 'Checked', 'off');
        
        set(shockTexts, 'Visible', 'on');
        
        if isempty(shockPoints.ZData)
            set(shockPoints, 'XData', T(:, 1),...
                'YData', T(:, 2));
        else
            set(shockPoints, 'XData', T(:, 1),...
                'YData', T(:, 2),...
                'ZData', T(:, 3));
        end
    end

    function selectShocksMenuCallback(~, ~)
        % dialog box
        validAnswer = false;
        while ~validAnswer
            dlgtitle = 'Shocks selection';
            prompt = {'Only display shocks #:'};
            dims = [1 35];
            definput = {''};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            if isempty(answer)
                return
            end
            
            visibleShocks = str2num(answer{1});
            if isempty(visibleShocks) || ~all(mod(visibleShocks, 1) == 0)...
                    || min(visibleShocks) < 1 || max(visibleShocks) > length(shockTexts)
                errorFig = errordlg('Incorrect input format', 'Error');
                waitfor(errorFig);
            else
                validAnswer = true;
            end
        end
        
        set(showAllMenu, 'Checked', 'off');
        set(selectShocksMenu, 'Checked', 'on');
        
        
        if isempty(shockPoints.ZData)
            set(shockPoints, 'XData', T(visibleShocks, 1),...
                'YData', T(visibleShocks, 2));
        else
            set(shockPoints, 'XData', T(visibleShocks, 1),...
                'YData', T(visibleShocks, 2),...
                'ZData', T(visibleShocks, 3));
        end
        
        set(shockTexts, 'Visible', 'off');
        set(shockTexts(visibleShocks), 'Visible', 'on');
        drawnow
    end

%% eigenvectors plot

pcNames = cell(1, Npc);
for k_pc = 1:Npc
    pcNames{k_pc} = sprintf('pc%u', k_pc);
end

if plotPC
    plotSpectrums(Freqs, W(:, 1:Npc)', spectrumFrequencyScale, spectrumScale, @(f) f/(2*QSpectrum),...
        ['Principal component basis', plotTitleSuffix], false, false, false, false, pcNames)
end

%% eigenvalues plot

if plotDistribution
    fig = figure('Name', ['Variance distribution', plotTitleSuffix]);
    ax = axes(fig);
    pie(ax, L(1:Npc) / sum(L));
    legend(ax, pcNames);
end


end

