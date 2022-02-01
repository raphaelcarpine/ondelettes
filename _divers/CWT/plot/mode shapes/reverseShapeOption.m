function shapePlot = reverseShapeOption(shapePlot0)
%REVERSESHAPEOPTION Summary of this function goes here
%   Detailed explanation goes here

    function shapePlot2(shape, varargin)
        fig = shapePlot0(shape, varargin{:});
        
        Rmenu = uimenu(fig, 'Text','REVERSE');
        function menuFct(~, ~)
            shapePlot2(-shape, varargin{:});
            delete(fig);
        end
        Rmenu.MenuSelectedFcn = @menuFct;
    end

shapePlot = @shapePlot2;

end

