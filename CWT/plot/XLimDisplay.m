classdef XLimDisplay < handle
    %XLIMDISPLAY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        XLim
        XLimRidges
        XLimVisible
        XLimRidgesVisible
        
        AxesXLim = gobjects(1, 0);
        updateLimAxesXLim = true(1, 0);
        AxesXLimRidges = gobjects(1, 0);
        updateLimAxesXLimRidges = true(1, 0);
        
        XminLines = gobjects(1, 0);
        XmaxLines = gobjects(1, 0);
        XminRidgesLines = gobjects(1, 0);
        XmaxRidgesLines = gobjects(1, 0);
    end
    
    methods
        function obj = XLimDisplay(XLim, XLimRidges, XLimVisible, XLimRidgesVisible, Axes, updateLimAxes)
            %XLIMDISPLAY Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 0 % test
                Axes = axes(figure);
                updateLimAxes = true;
                plot(Axes, 1:10, rand(3, 10));
                XLim = [2, 10];
                XLimRidges = [2, 7];
                XLimVisible = true;
                XLimRidgesVisible = true;
            end
            
            if nargin == 4
                Axes = gobjects(1, 0);
                updateLimAxes = true(1, 0);
            end
            
            obj.XLim = XLim;
            obj.XLimRidges = XLimRidges;
            obj.XLimVisible = XLimVisible;
            obj.XLimRidgesVisible = XLimRidgesVisible;
            
            obj.addAxes(Axes, updateLimAxes);
        end
        
        function addAxes(obj, Axes, updateLimAxes)
            obj.addAxesXLim(Axes, updateLimAxes);
            obj.addAxesXLimRidges(Axes, updateLimAxes);
        end
        
        function addAxesXLim(obj, Axes, updateLimAxes)
            Axes = Axes(:);
            updateLimAxes = updateLimAxes(:);
            for k_ax = 1:length(Axes)
                ax = Axes(k_ax);
                if ismember(ax, obj.AxesXLim)
                    continue
                end
                obj.AxesXLim(end+1) = ax;
                obj.updateLimAxesXLim(end+1) = updateLimAxes(k_ax);
                
                obj.XminLines(end+1) = xline(ax, obj.XLim(1), '-', 'Xmin',...
                    'LabelHorizontalAlignment', 'left', 'LabelOrientation', 'horizontal', 'Visible', obj.XLimVisible);
                obj.XmaxLines(end+1) = xline(ax, obj.XLim(2), '-', 'Xmax',...
                    'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Visible', obj.XLimVisible);
                obj.XminLines(end).Annotation.LegendInformation.IconDisplayStyle = 'off';
                obj.XmaxLines(end).Annotation.LegendInformation.IconDisplayStyle = 'off';
                
                drawnow
            end
        end
        
        function addAxesXLimRidges(obj, Axes, updateLimAxes)
            Axes = Axes(:);
            updateLimAxes = updateLimAxes(:);
            for k_ax = 1:length(Axes)
                ax = Axes(k_ax);
                if ismember(ax, obj.AxesXLimRidges)
                    continue
                end
                obj.AxesXLimRidges(end+1) = ax;
                obj.updateLimAxesXLimRidges(end+1) = updateLimAxes(k_ax);
                
                obj.XminRidgesLines(end+1) = xline(ax, obj.XLimRidges(1), '--', {'', 'Xmin ridges'},...
                    'LabelHorizontalAlignment', 'left', 'LabelOrientation', 'horizontal', 'Visible', obj.XLimRidgesVisible);
                obj.XmaxRidgesLines(end+1) = xline(ax, obj.XLimRidges(2), '--', {'', 'Xmax ridges'},...
                    'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Visible', obj.XLimRidgesVisible);
                obj.XminRidgesLines(end).Annotation.LegendInformation.IconDisplayStyle = 'off';
                obj.XmaxRidgesLines(end).Annotation.LegendInformation.IconDisplayStyle = 'off';
                
                drawnow
            end
        end
        
        function updateXLim(obj, XLim)
            obj.updateAxes();
            
            obj.XLim = XLim;
            for k_ax = 1:length(obj.AxesXLim)
                if ~obj.updateLimAxesXLim(k_ax)
                    continue
                end
                obj.XminLines(k_ax).Value = obj.XLim(1);
                obj.XmaxLines(k_ax).Value = obj.XLim(2);
            end
            
            drawnow
        end
        
        function updateXLimRidges(obj, XLimRidges)
            obj.updateAxes();
            
            obj.XLimRidges = XLimRidges;
            for k_ax = 1:length(obj.AxesXLimRidges)
                if ~obj.updateLimAxesXLimRidges(k_ax)
                    continue
                end
                obj.XminRidgesLines(k_ax).Value = obj.XLimRidges(1);
                obj.XmaxRidgesLines(k_ax).Value = obj.XLimRidges(2);
            end
            
            drawnow
        end
        
        function setVisible(obj, XLimVisible, XLimRidgesVisible)
            if nargin < 3
                XLimRidgesVisible = XLimVisible;
            end
            
            obj.updateAxes();
            
            obj.XLimVisible = XLimVisible;
            obj.XLimRidgesVisible = XLimRidgesVisible;
            set(obj.XminLines, 'Visible', obj.XLimVisible);
            set(obj.XmaxLines, 'Visible', obj.XLimVisible);
            set(obj.XminRidgesLines, 'Visible', obj.XLimRidgesVisible);
            set(obj.XmaxRidgesLines, 'Visible', obj.XLimRidgesVisible);
        end
    end
    
    methods (Access = private)
        function updateAxes(obj)
            obj.XminLines = obj.XminLines(isvalid(obj.AxesXLim));
            obj.XmaxLines = obj.XmaxLines(isvalid(obj.AxesXLim));
            obj.XminRidgesLines = obj.XminRidgesLines(isvalid(obj.AxesXLimRidges));
            obj.XmaxRidgesLines = obj.XmaxRidgesLines(isvalid(obj.AxesXLimRidges));
            obj.updateLimAxesXLim = obj.updateLimAxesXLim(isvalid(obj.AxesXLim));
            obj.updateLimAxesXLimRidges = obj.updateLimAxesXLimRidges(isvalid(obj.AxesXLimRidges));
            obj.AxesXLim = obj.AxesXLim(isvalid(obj.AxesXLim));
            obj.AxesXLimRidges = obj.AxesXLimRidges(isvalid(obj.AxesXLimRidges));
        end
    end
end

