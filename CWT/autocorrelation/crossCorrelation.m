function Rx = crossCorrelation(x)
%CROSSCORRELATION Summary of this function goes here
%   x(kDOF, kt)
Ndof = size(x, 1);
Nt = size(x, 2);
Rx = nan(Ndof, Ndof, Nt);

for i = 1:Ndof
    for j = i:Ndof
        Rxij = xcorr(x(i, :), x(j, :), 'biased');
        Rx(i, j, :) = Rxij(ceil(length(Rxij)/2):end);
        
        Rxji = flip(Rxij);
        Rx(j, i, :) = Rxji(ceil(length(Rxji)/2):end);
    end
end

end

