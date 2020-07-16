function [SVfftrx, SVfftvectrx] = svdFFT(Rx, Nsv)
%SVDCWT Summary of this function goes here
%   Rx : crossed correlation
%   Nsv : number of singular values
%   SVrx{k}(kf, kt) : k^th singular values of CWT(Rx) for each time and freq
%   SVvectrx{k}(:, kf, kt) : k^th singular vectors of CWT(Rx) for each time and freq

if nargin < 6
    Nsv = 1;
end

Ndof = size(Rx, 1);
Nt = size(Rx, 3);

%% CWT

FFTrx = nan(Ndof, Ndof, Nt);
for iddl = 1:Ndof
    for jddl = 1:Ndof
        FFTrx(iddl, jddl, :) = fft(reshape(Rx(iddl, jddl, :), [1, Nt]));
    end
end

%% SVD

SVfftrx = cell(1, Nsv);
SVfftvectrx = cell(1, Nsv);
for ksv = 1:Nsv %
    SVfftrx{ksv} = nan(1, Nt);
    SVfftvectrx{ksv} = nan(Ndof, Nt);
end


for kfreq = 1:Nt
    [U, S, ~] = svd(FFTrx(:, :, kfreq));
    for ksv = 1:Nsv
        SVfftrx{ksv}(kfreq) = S(ksv, ksv);
        SVfftvectrx{ksv}(:, kfreq) = U(:, ksv);
    end
    %TODO : phase
end

end

