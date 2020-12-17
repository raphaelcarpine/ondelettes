function [initFcn, updateFcn, closeFcn] = audioTimeOnAxes(Axes)
%AUDIOTIMEONAXES Summary of this function goes here
%   Detailed explanation goes here

timeLines = [];

    function initFcn0(t)
        for kax = 1:length(Axes)
            try
                hold(Axes(kax), 'on');
                timeLines(kax) = xline(Axes(kax), t, 'Color', '#4B0082');
            catch
                Axes = Axes([1:kax-1, kax+1:end]);
                timeLines = timeLines([1:kax-1, kax+1:end]);
            end
        end
        drawnow
    end

    function updateFcn0(t)
        for kax = 1:length(Axes)
            try
                set(timeLines(kax), 'value', t);
            catch
                Axes = Axes([1:kax-1, kax+1:end]);
                timeLines = timeLines([1:kax-1, kax+1:end]);
            end
        end
        drawnow
    end

    function closeFcn0()
        delete(timeLines);
        drawnow
    end


initFcn = @initFcn0;
updateFcn = @updateFcn0;
closeFcn = @closeFcn0;

end

