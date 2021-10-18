function [X, channelNames] = selectChannelsDirection(X, channelNames, direction)
%SELECTCHANNELSDIRECTION Summary of this function goes here
%   Detailed explanation goes here

I = cellfun(@(ch) ch(end) == direction, channelNames);

X = X(I, :);
channelNames = channelNames(I);

end

