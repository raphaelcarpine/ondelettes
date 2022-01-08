clear all
close all




%% système

% le système est le suivant
%
%             ______________                               ______________
%   \|   ____| ressort k(1) |____     ____________    ____| ressort k(2) |____     ____________
%   \|__|    |______________|    |___| masse m(1) |__|    |______________|    |___| masse m(2) |
%   \|  |   __________________   |   |____________|  |   __________________   |   |____________|
%   \|  |__| amortisseur c(1) |__|                   |__| amortisseur c(2) |__|
%          |__________________|                         |__________________|
%

m = [1, 1];
k = 100 * [1, 1];
c = [0.5, 0];
k = [k, 0];
c = [c, 0];

M = diag(m);
K = zeros(2);
for ik = 1:3
    if ik - 1 > 0
        K(ik-1, ik-1) = K(ik-1, ik-1) + k(ik);
    end
    if ik <= 2
        K(ik, ik) = K(ik, ik) + k(ik);
    end
    if ik - 1 > 0 && ik <= 2
        K(ik, ik-1) = K(ik, ik-1) - k(ik);
        K(ik-1, ik) = K(ik-1, ik) - k(ik);
    end
end
C = zeros(2);
for ic = 1:3
    if ic - 1 > 0
        C(ic-1, ic-1) = C(ic-1, ic-1) + c(ic);
    end
    if ic <= 2
        C(ic, ic) = C(ic, ic) + c(ic);
    end
    if ic - 1 > 0 && ic <= 2
        C(ic, ic-1) = C(ic, ic-1) - c(ic);
        C(ic-1, ic) = C(ic-1, ic) - c(ic);
    end
end

syst = systLin(M, K, C);

poles = syst.complexModes();
freqs_propres = abs(poles)/(2*pi)
amorts = -100*real(poles)./abs(poles)
mu = -real(poles);
freqs_amort = abs(imag(poles))/(2*pi);
[Mb, Kb, Cb] = syst.modalDamping();
% Kb^(-1/2)*Cb

T = 100;
dt = 0.05;
fe = 1/dt;

nt = round(T/dt);
t = dt * (0:nt-1);
T = nt*dt;

mu_tilde = nan(2, 0);

ti = nan(1, 2);

N_sim = 100;
[initWaitBar, updateWaitBar, closeWaitBar] = getWaitBar(N_sim);
initWaitBar();
for k_sim = 1:N_sim
    % excitation
    f = randn(1, nt);
    
    % reponse
    ddlF = 1;
    x = syst.response(f, dt, ddlF);
    
    % autocorrelation
    Rx = xcorr(x(ddlF, :), 'unbiased');
    Rx = Rx(ceil(length(Rx)/2):end);
    
    for k_mode = 1:2
        % ondelette
        ct = 3;
        cf = 5;
        MotherWavelet = 'cauchy';
        [Qmin, Qmax, Qz] = getBoundsQ2(freqs_amort(k_mode), diff(freqs_amort),...
            amorts(k_mode), 0, 0, ct, ct, cf, MotherWavelet);
        Q = Qmin + 1;
        fmin = freqs_amort(k_mode) - 0.5;
        fmax = freqs_amort(k_mode) + 0.5;
        
        Topt = 1.4326/mu(k_mode);
        [~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
        ti(k_mode) = DeltaT(freqs_amort(k_mode));
        
        Rx = Rx(t <= Topt + 2*ti(k_mode) + 1);
        tRx = t(t <= Topt + 2*ti(k_mode) + 1);
        
%         % wavelet menu
%         figure;
%         plt = plot(tRx, Rx);
%         WaveletMenu('WaveletPlot', plt, 'Q', Q, 'fmin', fmin, 'fmax', fmax,...
%             'MotherWavelet', MotherWavelet, 'XLimRidge', [0, ti + Topt]);
        
        % ridge
        ridge = RidgeExtract(tRx, Rx, Q, fmin, fmax, 100, 'NbMaxRidges', 1,...
            'MotherWavelet', MotherWavelet, 'XLimRidge', [0, ti(k_mode) + Topt]);
        if length(ridge.freq) ~= 1
            error(' ');
        end
        
        % regression
        logAmpl = log(abs(ridge.val{1}));
        ridgeTime = ridge.time{1};
        reg_coeffs = [ones(size(ridgeTime)); ridgeTime].' \ logAmpl.';
        mu_tilde(k_mode, k_sim) = -reg_coeffs(2);
        
%         % plot
%         figure;
%         plot(ridgeTime, ridge.freq{1});
%         
%         figure;
%         plot(ridgeTime, logAmpl);
%         hold on
%         plot(ridgeTime, [ones(size(ridgeTime)); ridgeTime].' * reg_coeffs);
%         xline(ti + Topt, '--');
        
    end
    updateWaitBar(k_sim);
end
closeWaitBar();

%%

for k_mode = 1:2
    mu_tilde_mean = mean(mu_tilde(k_mode, :));
    mu_tilde_std = std(mu_tilde(k_mode, :));
    alpha_p = sqrt(4/(mu(k_mode)*T) + 2/(freqs_amort(k_mode)*T));
    Meff = mu_tilde_std/(mu(k_mode) * alpha_p * exp(mu(k_mode)*ti(k_mode)));
    fprintf('mode %d (mu=%.4f):\nmu*T=%.1f, alpha=%.3f\nmean: %.4f, std dev: %.4f\nMeff: %.2f\n\n',...
        [k_mode, mu(k_mode), mu(k_mode)*T, alpha_p, mu_tilde_mean, mu_tilde_std, Meff]);
end











