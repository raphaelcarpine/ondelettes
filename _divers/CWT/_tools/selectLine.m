function selectLine(ax)
%SELECTLINE Summary of this function goes here
%   Detailed explanation goes here
if nargin == 0
    ax = gca;
end

Lines = findobj(ax, 'Type', 'line');
LineWidths = [Lines.LineWidth];

    function lineCallback(hObject, eventData)
        0;
        for k = 1:length(Lines)
            if Lines(k) ~= hObject
                if eventData.Button == 1
                    Lines(k).LineWidth = LineWidths(k);
                end
            elseif Lines(k).LineWidth == LineWidths(k)
                Lines(k).LineWidth = 3*LineWidths(k);
            else
                Lines(k).LineWidth = LineWidths(k);
            end
        end
    end

set(Lines, 'buttonDownFcn', @lineCallback);

end

