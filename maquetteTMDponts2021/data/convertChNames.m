function [X, chNames, captNumber, captDir] = convertChNames(X, chNames)
%CONVERTCHNAMES Summary of this function goes here
%   Detailed explanation goes here

for kch = 1:length(chNames)
    switch str2double(chNames{kch}(end))
        case 1
            chNames{kch}(end) = 'x';
        case 2
            chNames{kch}(end) = 'y';
            X(kch, :) = -X(kch, :);
        case 3
            chNames{kch}(end) = 'z';
            X(kch, :) = -X(kch, :);
        otherwise
            error('');
    end
end

[chNames, I] = sort(chNames);
X = X(I, :);


captNumber = [];
captDir = {};
for kch = 1:length(chNames)
    chName = chNames{kch};
    if isempty(captNumber) || str2double(chName(1:5)) ~= captNumber(end)
        captNumber(end+1) = str2double(chName(1:5));
        captDir{end+1} = chName(end);
    else
        captDir{end}(end+1) = chName(end);
    end
end

end

