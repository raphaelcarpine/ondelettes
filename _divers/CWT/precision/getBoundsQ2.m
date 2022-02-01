function [Qmin, Qmax, Qz] = getBoundsQ2(f, Df, zeta, TLim, TLimRidge, ct, ctRidge, cf, MotherWavelet)
%GETBOUNDSQ Le, Argoul 2004
%   f = (f1+f2)/2
%   Df = f2 - f1
%   Dt = 1/lambda
%   T = total time

%% Qmin
Qmin = cf * f/(2*Df);

%% Qmax

    function t = remainingTime(Q)
        [~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
        Deltat = DeltaT(f);
        t1 = TLim + Deltat * [1 -1].*ct;
        t2 = TLimRidge + Deltat * [1 -1].*ctRidge;
        t = min(t1(2), t2(2)) - max(t1(1), t2(1));
    end

Qmax = fminsearch(@(Q) abs(remainingTime(Q)), Qmin);

% Qmax = 2*pi*f*T / (4*ct*1/2); %mu_psi \simeq 1/2

%% Qz
%Qz = 1/sqrt(2) / zeta / ct; %autre def
Qz = 1/sqrt(2) / zeta;
end

