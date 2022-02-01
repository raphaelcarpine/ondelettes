function [WvltOut, IndWvltOut] = raindrop(WvltIn)

WvltOut = abs(WvltIn); % On ne travaille qu'en module

WvltOut_inf = [WvltOut(1, :) + 1; WvltOut(1:end-1, :)]; % pas de maximum locaux aux bords
WvltOut_sup = [WvltOut(2:end, :); WvltOut(end, :) + 1]; % pas de maximum locaux aux bords

IndWvltOut = WvltOut >= WvltOut_inf & WvltOut > WvltOut_sup;

WvltOut(~IndWvltOut) = nan;
end
