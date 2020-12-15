function [freqs, SVfftrx, SVfftvectrx] = svdFFT(tRx, Rx, Nsv, varargin)
%SVDCWT Summary of this function goes here
%   Rx : crossed correlation
%   Nsv : number of singular values
%   SVfftrx{k}(kf) : k^th singular values of CWT(Rx) for each time and freq
%   SVvectrx{k}(:, kf) : k^th singular vectors of CWT(Rx) for each time and freq

Ndof = size(Rx, 1);
Nt = size(Rx, 3);

%% CWT

Nf = floor(Nt/2);

FFTrx = nan(Ndof, Ndof, Nf);
for iddl = 1:Ndof
    for jddl = 1:Ndof
        [freqs, FFTrx(iddl, jddl, :)] = fourierTransform(tRx, reshape(Rx(iddl, jddl, :), [1, Nt]), varargin{:});
    end
end

%% SVD

SVfftrx = cell(1, Nsv);
SVfftvectrx = cell(1, Nsv);
for ksv = 1:Nsv %
    SVfftrx{ksv} = nan(1, Nf);
    SVfftvectrx{ksv} = nan(Ndof, Nf);
end


for kfreq = 1:Nf
    [U, S, ~] = svd(FFTrx(:, :, kfreq));
    for ksv = 1:Nsv
        SVfftrx{ksv}(kfreq) = S(ksv, ksv);
        SVfftvectrx{ksv}(:, kfreq) = U(:, ksv);
    end
    %TODO : phase
end

end

