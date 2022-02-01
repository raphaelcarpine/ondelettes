classdef AudioTimeOnAxes < handle
    %AUDIOTIMEONAXES2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        timeLines = gobjects(1, 0);
        Axes = gobjects(1, 0);
        t
    end
    
    methods
        function obj = AudioTimeOnAxes(t, Axes)
            %AUDIOTIMEONAXES2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.t = t;
            
            if nargin > 1
                obj.addAxes(Axes);
            end
            drawnow
        end
        
        function addAxes(obj, Axes)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            for kax = 1:length(Axes)
                try
                    hold(Axes(kax), 'on');
                    obj.timeLines(end+1) = xline(Axes(kax), obj.t, 'Color', '#4B0082');
                    obj.timeLines(end).Annotation.LegendInformation.IconDisplayStyle = 'off';
                    obj.Axes(end+1) = Axes(kax);
                catch
                end
            end
        end
        
        function updateT(obj, t)
            obj.t = t;
            kax = 1;
            while kax <= length(obj.Axes)
                try
                    set(obj.timeLines(kax), 'value', obj.t);
                    kax = kax+1;
                catch
                    obj.Axes = obj.Axes([1:kax-1, kax+1:end]);
                    obj.timeLines = obj.timeLines([1:kax-1, kax+1:end]);
                end
            end
            drawnow
        end
        
        function close(obj)
            delete(obj.timeLines);
            drawnow
        end
    end
end

