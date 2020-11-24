function [initFcn, updateFcn, closeFcn] = audioTimeOnAxes(Axes)
%AUDIOTIMEONAXES Summary of this function goes here
%   Detailed explanation goes here

Nax = length(Axes);

timeLines = [];

    function initFcn0(t)
        for kax = 1:Nax
            hold(Axes(kax), 'on');
            timeLines(kax) = xline(Axes(kax), t, 'Color', '#4B0082');
        end
        drawnow
    end

    function updateFcn0(t)
        for kax = 1:Nax
            set(timeLines(kax), 'value', t);
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

