classdef PlotRDTtimes < handle
    %PLOTRDTTIMES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        t
        X
        removeMean
        referenceChannel
        ax
        linePlots
        scatterPlot
        thresholdLine
        lgdTitle
    end
    
    methods
        function obj = PlotRDTtimes(t, X, removeMean, referenceChannel, signalChannels)
            %PLOTRDTTIMES Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                removeMean = true;
            end
            if nargin < 4
                referenceChannel = 1;
            end
            if nargin < 5
                signalChannels = 1:size(X, 1);
            end
            
            obj.t = t;
            obj.X = X;
            obj.removeMean = removeMean;
            obj.referenceChannel = referenceChannel;
            
            fig = figure('Name', 'RDT times', 'NumberTitle','off');
            obj.ax = axes(fig);
            hold(obj.ax, 'on');
            
            yline(obj.ax, 0, '--');
            obj.linePlots = plot(obj.ax, t, X - obj.removeMean * mean(X, 2), 'LineStyle', ':');
            set(obj.linePlots(obj.referenceChannel), 'LineStyle', '-');
            obj.thresholdLine = yline(obj.ax, 0, '-r', 'threshold', 'LabelHorizontalAlignment', 'left',...
                'Visible', false);
            obj.scatterPlot = scatter(obj.ax, [], [], '+b', 'LineWidth', 2);
            
            % legend
            legendNames = cell(1, size(X, 1));
            for k_ch = 1:size(X, 1)
                legendNames{k_ch} = sprintf('channel %u', signalChannels(k_ch));
            end
            lgd = legend(obj.linePlots, legendNames);
            obj.lgdTitle = title(lgd, 'Test');
        end
        
        function update(obj, removeMean, referenceChannel, threshold, thresholdString, KTrdt)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if ~isvalid(obj.ax)
                obj = PlotRDTtimes(obj.t, obj.X, obj.removeMean, obj.referenceChannel);
            end
            
            % update
            if removeMean ~= obj.removeMean
                obj.removeMean = removeMean;
                for k = 1:length(obj.linePlots)
                    set(obj.linePlots(k), 'YData', obj.X(k, :) - obj.removeMean * mean(obj.X(k, :)));
                end
            end
            
            if referenceChannel ~= obj.referenceChannel
                set(obj.linePlots(obj.referenceChannel), 'LineStyle', ':');
                obj.referenceChannel = referenceChannel;
                set(obj.linePlots(obj.referenceChannel), 'LineStyle', '-');
            end
            
            set(obj.thresholdLine, 'Value', threshold, 'Visible', true);
            set(obj.thresholdLine, 'Label', ['threshold = ', thresholdString]);
            set(obj.scatterPlot, 'XData', obj.t(KTrdt), 'YData', threshold*ones(size(KTrdt)));
            set(obj.lgdTitle, 'String', sprintf('Nb of detections: %u', length(KTrdt)));
            drawnow
        end
    end
end

