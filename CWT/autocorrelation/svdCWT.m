function [SVrx, SVvectrx] = svdCWT(t, Rx, WvltFreqs, Q, Nsv)
%SVDCWT Summary of this function goes here
%   Rx : crossed correlation
%   Nsv : number of singular values
%   SVrx{k}(kf, kt) : k^th singular values of CWT(Rx) for each time and freq
%   SVvectrx{k}(:, kf, kt) : k^th singular vectors of CWT(Rx) for each time and freq

timeWaitBarSVD = 5;
timeWaitBarCWT = 5;

if nargin < 6
    Nsv = 1;
end

Ndof = size(Rx, 1);
Nt = size(Rx, 3);
NbFreq = length(WvltFreqs);

%% CWT

CWTrx = nan(Ndof, Ndof, NbFreq, Nt);

tic;
for iddl = 1:Ndof
    for jddl = 1:Ndof
        r = (Ndof*(iddl-1) + jddl)/Ndof^2;
        if ~isnan(timeWaitBarCWT) && toc > timeWaitBarCWT
            w = waitbar(r, [num2str(round(100*r)), '%'], 'Name', 'Computing CWT');
            timeWaitBarCWT = nan;
        elseif isnan(timeWaitBarCWT) && isvalid(w)
            waitbar(r, w, [num2str(round(100*r)), '%']);
        end
        
        CWTrx(iddl, jddl, :, :) = WvltComp(t, reshape(Rx(iddl, jddl, :), [1, Nt]), WvltFreqs, Q);
    end
end

if isnan(timeWaitBarCWT) && isvalid(w)
    waitbar(1, w, '100%');
    close(w);
end

%% SVD

SVrx = cell(1, Nsv);
SVvectrx = cell(1, Nsv);
for ksv = 1:Nsv %
    SVrx{ksv} = nan(NbFreq, Nt);
    SVvectrx{ksv} = nan(Ndof, NbFreq, Nt);
end


tic;
for kfreq = 1:NbFreq
    r = (kfreq-1)/NbFreq;
    if ~isnan(timeWaitBarSVD) && toc > timeWaitBarSVD
        w = waitbar(r, [num2str(round(100*r)), '%'], 'Name', 'Computing SVD');
        timeWaitBarSVD = nan;
    elseif isnan(timeWaitBarSVD) && isvalid(w)
        waitbar(r, w, [num2str(round(100*r)), '%']);
    end
    
    for kt = 1:Nt
        [U, S, ~] = svd(CWTrx(:, :, kfreq, kt));
        for ksv = 1:Nsv
            SVrx{ksv}(kfreq, kt) = S(ksv, ksv);
            SVvectrx{ksv}(:, kfreq, kt) = U(:, ksv);
        end
        %TODO : phase
    end
end

if isnan(timeWaitBarSVD) && isvalid(w)
    waitbar(1, w, '100%');
    close(w);
end

end

