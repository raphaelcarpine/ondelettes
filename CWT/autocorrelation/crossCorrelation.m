function Rx = crossCorrelation(x, maxLag)
%CROSSCORRELATION Summary of this function goes here
%   x(kDOF, kt)
Ndof = size(x, 1);
Rx = nan(Ndof, Ndof, maxLag+1);

for i = 1:Ndof
    for j = i:Ndof
        Rxij = xcorr(x(i, :), x(j, :), maxLag, 'biased');
        Rx(i, j, :) = Rxij(ceil(length(Rxij)/2):end);
        
        Rxji = flip(Rxij);
        Rx(j, i, :) = Rxji(ceil(length(Rxji)/2):end);
    end
end

end

