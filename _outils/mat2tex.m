function s = mat2tex(mat, formatSpec)
%MAT2TEX Summary of this function goes here
%   Detailed explanation goes here

s = '\begin{tabular}{';
s = [s, repmat('c', 1, size(mat, 2)), '}', newline];

for i = 1:size(mat, 1)
    for j = 1:size(mat, 2)
        s = [s, sprintf(formatSpec, mat(i, j))];
        if j < size(mat, 2)
            s = [s, ' & '];
        elseif i < size(mat, 1)
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

