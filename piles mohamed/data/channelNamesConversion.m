function [X, channelNames] = channelNamesConversion(X, channelNames)
%CHANNELNAMESCONVERSION Summary of this function goes here
%   Detailed explanation goes here

X = X.';

for kch = 1:length(channelNames)
    chName = channelNames{kch};
    
    % accelero
    switch str2double(chName(1:5))
        case 29277
            chAcc = 1;
        case 40202
            chAcc = 2;
        case 40197
            chAcc = 3;
        case 40201
            chAcc = 4;
        case 40200
            chAcc = 5;
        otherwise
            error('');
    end
    
    % axis
    if chAcc == 1
        switch str2double(chName(end))
            case 1
                chAx = 'y';
            case 2
                chAx = 'x';
            case 3
                chAx = 'z';
            otherwise
                error('');
        end
    else
        switch str2double(chName(end))
            case 1
                chAx = 'z';
            case 2
                chAx = 'y';
            case 3
                chAx = 'x';
            otherwise
                error('');
        end
    end
    
    channelNames{kch} = ['acc', num2str(chAcc), chAx];
    
end

[channelNames, I] = sort(channelNames);
X = X(I, :);

X(1, :) = -X(1, :); % acc1x opposé

X = X.';

end

