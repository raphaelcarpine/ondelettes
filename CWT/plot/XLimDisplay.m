classdef XLimDisplay < handle
    %XLIMDISPLAY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        XLim
        XLimRidges
        XLimVisible
        XLimRidgesVisible
        
        Axes = gobjects(1, 0);
        updateLimAxes = true(1, 0);
        
        XminLines = gobjects(1, 0);
        XminTxts = gobjects(1, 0);
        XmaxLines = gobjects(1, 0);
        XmaxTxts = gobjects(1, 0);
        XminRidgesLines = gobjects(1, 0);
        XminRidgesTxts = gobjects(1, 0);
        XmaxRidgesLines = gobjects(1, 0);
        XmaxRidgesTxts = gobjects(1, 0);
    end
    
    methods
        function obj = XLimDisplay(XLim, XLimRidges, XLimVisible, XLimRidgesVisible, Axes, updateLimAxes)
            %XLIMDISPLAY Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 0 % test
                Axes = axes;
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
            Axes = Axes(:);
            updateLimAxes = updateLimAxes(:);
            for k_ax = 1:length(Axes)
                ax = Axes(k_ax);
                obj.Axes(end+1) = ax;
                obj.updateLimAxes(end+1) = updateLimAxes(k_ax);
                
                obj.XminLines(end+1) = xline(ax, obj.XLim(1), 'Visible', obj.XLimVisible);
                obj.XmaxLines(end+1) = xline(ax, obj.XLim(2), 'Visible', obj.XLimVisible);
                obj.XminRidgesLines(end+1) = xline(ax, obj.XLimRidges(1), 'Visible', obj.XLimRidgesVisible, 'LineStyle', '--');
                obj.XmaxRidgesLines(end+1) = xline(ax, obj.XLimRidges(2), 'Visible', obj.XLimRidgesVisible, 'LineStyle', '--');
                
                text1 = text(ax, obj.XLim(1), ax.YLim(2), 'Xmin ', 'Visible', obj.XLimVisible,...
                    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Clipping', 'on');
                text2 = text(ax, obj.XLim(2), ax.YLim(2), ' Xmax', 'Visible', obj.XLimVisible,...
                    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Clipping', 'on');
                text3 = text(ax, obj.XLimRidges(1), ax.YLim(2), {'', 'Xmin ridges '}, 'Visible', obj.XLimRidgesVisible,...
                    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Clipping', 'on');
                text4 = text(ax, obj.XLimRidges(2), ax.YLim(2), {'', ' Xmax ridges'}, 'Visible', obj.XLimRidgesVisible,...
                    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Clipping', 'on');
                obj.XminTxts(end+1) = text1;
                obj.XmaxTxts(end+1) = text2;
                obj.XminRidgesTxts(end+1) = text3;
                obj.XmaxRidgesTxts(end+1) = text4;
                
                % vertical alignment with height of axes
                addlistener(ax, 'YLim', 'PostSet',...
                    @(~,~) set(text1, 'Position', get(text1, 'Position')*[1;0;0] * [1 0] + get(ax, 'YLim').*[0 1]));
                addlistener(ax, 'YLim', 'PostSet',...
                    @(~,~) set(text2, 'Position', get(text2, 'Position')*[1;0;0] * [1 0] + get(ax, 'YLim').*[0 1]));
                addlistener(ax, 'YLim', 'PostSet',...
                    @(~,~) set(text3, 'Position', get(text3, 'Position')*[1;0;0] * [1 0] + get(ax, 'YLim').*[0 1]));
                addlistener(ax, 'YLim', 'PostSet',...
                    @(~,~) set(text4, 'Position', get(text4, 'Position')*[1;0;0] * [1 0] + get(ax, 'YLim').*[0 1]));
                
                drawnow
            end
        end
        
        function updateXLim(obj, XLim)
            obj.updateAxes();
            
            obj.XLim = XLim;
            for k_ax = 1:length(obj.Axes)
                if ~obj.updateLimAxes(k_ax)
                    continue
                end
                obj.XminLines(k_ax).Value = obj.XLim(1);
                obj.XminTxts(k_ax).Position(1) = obj.XLim(1);
                obj.XmaxLines(k_ax).Value = obj.XLim(2);
                obj.XmaxTxts(k_ax).Position(1) = obj.XLim(2);
            end
            
            drawnow
        end
        
        function updateXLimRidges(obj, XLimRidges)
            obj.updateAxes();
            
            obj.XLimRidges = XLimRidges;
            for k_ax = 1:length(obj.Axes)
                if ~obj.updateLimAxes(k_ax)
                    continue
                end
                obj.XminRidgesLines(k_ax).Value = obj.XLimRidges(1);
                obj.XminRidgesTxts(k_ax).Position(1) = obj.XLimRidges(1);
                obj.XmaxRidgesLines(k_ax).Value = obj.XLimRidges(2);
                obj.XmaxRidgesTxts(k_ax).Position(1) = obj.XLimRidges(2);
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
            set(obj.XminTxts, 'Visible', obj.XLimVisible);
            set(obj.XmaxLines, 'Visible', obj.XLimVisible);
            set(obj.XmaxTxts, 'Visible', obj.XLimVisible);
            set(obj.XminRidgesLines, 'Visible', obj.XLimRidgesVisible);
            set(obj.XminRidgesTxts, 'Visible', obj.XLimRidgesVisible);
            set(obj.XmaxRidgesLines, 'Visible', obj.XLimRidgesVisible);
            set(obj.XmaxRidgesTxts, 'Visible', obj.XLimRidgesVisible);
        end
    end
    
    methods (Access = private)
        function updateAxes(obj)
            obj.XminLines = obj.XminLines(isvalid(obj.Axes));
            obj.XminTxts = obj.XminTxts(isvalid(obj.Axes));
            obj.XmaxLines = obj.XmaxLines(isvalid(obj.Axes));
            obj.XmaxTxts = obj.XmaxTxts(isvalid(obj.Axes));
            obj.XminRidgesLines = obj.XminRidgesLines(isvalid(obj.Axes));
            obj.XminRidgesTxts = obj.XminRidgesTxts(isvalid(obj.Axes));
            obj.XmaxRidgesLines = obj.XmaxRidgesLines(isvalid(obj.Axes));
            obj.XmaxRidgesTxts = obj.XmaxRidgesTxts(isvalid(obj.Axes));
            obj.updateLimAxes = obj.updateLimAxes(isvalid(obj.Axes));
            obj.Axes = obj.Axes(isvalid(obj.Axes));
        end
    end
end

