clear all
close all




%% système

% le système est le suivant
%
%             ______________
%   \|   ____| ressort k(1) |____     ____________
%   \|__|    |______________|    |___| masse m(1) |
%   \|  |   __________________   |   |____________|
%   \|  |__| amortisseur c(1) |__|
%          |__________________|
%

m = 1;
k = 100;
c = 0.5;

syst = systLin(m, k, c);

poles = syst.complexModes();
freqs_propres = abs(poles)/(2*pi)
amorts = -100*real(poles)./abs(poles)
mu = -real(poles);
freqs_amort = abs(imag(poles))/(2*pi);
[Mb, Kb, Cb] = syst.modalDamping();
% Kb^(-1/2)*Cb

T = 10000;
dt = 0.05;
fe = 1/dt;

nt = round(T/dt);
t = dt * (0:nt-1);
T = nt*dt;

mu_tilde = nan(1, 0);

ti = nan;

N_sim = 1000;
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
    
    % ondelette
    ct = 3;
    cf = 5;
    MotherWavelet = 'cauchy';
    Q = 2;
    fmin = 1;
    fmax = 2;
    
    Topt = 1.4326/mu;
    [~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
    ti = DeltaT(freqs_amort);
    
    Rx = Rx(t <= Topt + 2*ti + 1);
    tRx = t(t <= Topt + 2*ti + 1);
    
    %         % wavelet menu
    %         figure;
    %         plt = plot(tRx, Rx);
    %         WaveletMenu('WaveletPlot', plt, 'Q', Q, 'fmin', fmin, 'fmax', fmax,...
    %             'MotherWavelet', MotherWavelet, 'XLimRidge', [0, ti + Topt]);
    
    % ridge
    ridge = RidgeExtract(tRx, Rx, Q, fmin, fmax, 100, 'NbMaxRidges', 1,...
        'MotherWavelet', MotherWavelet, 'XLimRidge', [0, ti + Topt]);
    if length(ridge.freq) ~= 1
        error(' ');
    end
    
    % regression
    logAmpl = log(abs(ridge.val{1}));
    ridgeTime = ridge.time{1};
    reg_coeffs = [ones(size(ridgeTime)); ridgeTime].' \ logAmpl.';
    mu_tilde(k_sim) = -reg_coeffs(2);
    
    %         % plot
    %         figure;
    %         plot(ridgeTime, ridge.freq{1});
    %
    %         figure;
    %         plot(ridgeTime, logAmpl);
    %         hold on
    %         plot(ridgeTime, [ones(size(ridgeTime)); ridgeTime].' * reg_coeffs);
    %         xline(ti + Topt, '--');
    
    updateWaitBar(k_sim);
end
closeWaitBar();

%%

mu_tilde_mean = mean(mu_tilde(:));
mu_tilde_std = std(mu_tilde(:));
alpha_p = sqrt(4/(mu*T) + 2/(freqs_amort*T));
Meff = mu_tilde_std/(mu * alpha_p * exp(mu*ti));
fprintf('mode %d (mu=%.4f):\nmu*T=%.1f, alpha=%.3f\nmean: %.4f, std dev: %.4f\nMeff: %.2f\n\n',...
    [1, mu, mu*T, alpha_p, mu_tilde_mean, mu_tilde_std, Meff]);












