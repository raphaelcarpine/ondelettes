function Rx = crossCorrelation(x, maxLag)
%CROSSCORRELATION Summary of this function goes here
%   x(kDOF, kt)

timeWaitBar = 5;

Ndof = size(x, 1);
Rx = nan(Ndof, Ndof, maxLag+1);

tic;
for i = 1:Ndof
    for j = i:Ndof
        r = (i*(Ndof-1)+j-1)/Ndof^2;
        if ~isnan(timeWaitBar) && toc > timeWaitBar
            w = waitbar(r, [num2str(round(100*r)), '%'], 'Name', 'computing cross-corr');
            timeWaitBar = nan;
        elseif isnan(timeWaitBar) && isvalid(w)
            waitbar(r, w, [num2str(round(100*r)), '%']);
        end
        
        Rxij = xcorr(x(i, :), x(j, :), maxLag, 'biased');
        Rx(i, j, :) = Rxij(ceil(length(Rxij)/2):end);
        
        Rxji = flip(Rxij);
        Rx(j, i, :) = Rxji(ceil(length(Rxji)/2):end);
    end
end

if isnan(timeWaitBar) && isvalid(w)
    waitbar(1, w, '100%');
    close(w);
end

end

