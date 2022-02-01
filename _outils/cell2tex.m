function s = cell2tex(C, formatSpec)
%MAT2TEX Summary of this function goes here
%   Detailed explanation goes here

s = '\begin{tabular}{';
s = [s, repmat('c', 1, size(C, 2)), '}', newline];

for i = 1:size(C, 1)
    for j = 1:size(C, 2)
        if ischar(C{i, j})
            s = [s, C{i, j}];
        else
            s = [s, sprintf(formatSpec, C{i, j})];
        end
        if j < size(C, 2)
            s = [s, ' & '];
        elseif i < size(C, 1)
            s = [s, ' \\', newline];
        else
            s = [s, newline];
        end
    end
end
s = [s, '\end{tabular}', newline];

if nargout == 0
    fprintf('%s', s);
end

end

