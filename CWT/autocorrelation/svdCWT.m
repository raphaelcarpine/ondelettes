function [SVrx, SVvectrx] = svdCWT(t, Rx, WvltFreqs, Q, Nsv, varargin)
%SVDCWT Summary of this function goes here
%   Rx : crossed correlation
%   Nsv : number of singular values
%   SVrx{k}(kf, kt) : k^th singular values of CWT(Rx) for each time and freq
%   SVvectrx{k}(:, kf, kt) : k^th singular vectors of CWT(Rx) for each time and freq

timeWaitBarSVD = 5;
timeWaitBarCWT = 5;

if nargin < 5
    Nsv = 1;
end

Ndof = size(Rx, 1);
Nt = size(Rx, 3);
NbFreq = length(WvltFreqs);

%% CWT

CWTrx = nan(Ndof, Ndof, NbFreq, Nt);

[initWaitBar, updateWaitBar, closeWaitBar] =...
    getWaitBar(Ndof^2, 'windowTitle', 'Computing CWTs', 'displayTime', 0);
initWaitBar();
for iddl = 1:Ndof
    for jddl = 1:Ndof
        updateWaitBar();
        
        CWTrx(iddl, jddl, :, :) = WvltComp(t, reshape(Rx(iddl, jddl, :), [1, Nt]),...
            WvltFreqs, Q, varargin{:});
    end
end
closeWaitBar();

%% SVD

SVrx = cell(1, Nsv);
SVvectrx = cell(1, Nsv);
for ksv = 1:Nsv %
    SVrx{ksv} = nan(NbFreq, Nt);
    SVvectrx{ksv} = nan(Ndof, NbFreq, Nt);
end

[initWaitBar, updateWaitBar, closeWaitBar] =...
    getWaitBar(NbFreq, 'windowTitle', 'Computing SVDs', 'displayTime', 0);
initWaitBar();
for kfreq = 1:NbFreq
    updateWaitBar();
    
    for kt = 1:Nt
        [U, S, ~] = svd(CWTrx(:, :, kfreq, kt));
        for ksv = 1:Nsv
            SVrx{ksv}(kfreq, kt) = S(ksv, ksv);
            SVvectrx{ksv}(:, kfreq, kt) = U(:, ksv);
        end
        %TODO : phase
    end
end

closeWaitBar();

end

