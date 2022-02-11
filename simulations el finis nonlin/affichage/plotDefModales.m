function plotDefModales(dispStresses, dispFreqs, plotShapes, L, E, J, mu, N, Mr, Cr, Kr, fleche_pont, sigma0, sigma_min, sigma_max, delta_sigma_1t, delta_sigma_1t_12m)

% precontrainte etc.
if dispStresses
    fprintf('max deflection static: %.1f mm (linear)\n', 1000*fleche_pont);
    fprintf('pre-stress: %.1f MPa\n', sigma0*1e-6);
    fprintf('min. compression: %.2f MPa\n', sigma_min*1e-6);
    fprintf('max. compression: %.1f MPa\n', sigma_max*1e-6);
    fprintf('min. compression loss 1t vehicle: %.2f MPa\n', delta_sigma_1t*1e-6);
    fprintf('min. compression loss 2t vehicle: %.2f MPa\n', 2*delta_sigma_1t*1e-6);
    fprintf('min. compression loss 20t vehicle: %.2f MPa\n', 20*delta_sigma_1t_12m*1e-6);
    fprintf('min. compression loss 30t vehicle: %.2f MPa\n', 30*delta_sigma_1t_12m*1e-6);
    fprintf('min. compression loss 40t vehicle: %.2f MPa\n', 40*delta_sigma_1t_12m*1e-6);
end

% frequences
if dispFreqs
    % frequences therorique, problème continu
    freqsTh = pi/(2*L^2) * sqrt(E*J/mu) * (1:N-2).^2;
    
    % frequence propre problème discrétisé
    freqs = eig(Mr\Kr);
    freqs = sqrt(freqs)/(2*pi);
    freqs = sort(freqs);
    
    % poles
    poles = eig([zeros(size(Mr)), eye(size(Mr)); -Mr\Kr, -Mr\Cr]);
    poles = poles(imag(poles) ~= 0);
    [~, Ipoles] = sort(imag(poles));
    poles = poles(Ipoles(end/2+1:end));
    amorts = -real(poles)./abs(poles);
    amorts = [amorts.', inf*ones(1, length(freqs)-length(amorts))];
    
    % affichage
    for kfreq = 1:5
        fprintf('mode %d: f=%.3fHz (f_th=%.3fHz, %.2f%% error), z=%.2f%%\n',...
            [kfreq, freqs(kfreq), freqsTh(kfreq), (freqs(kfreq)-freqsTh(kfreq))/freqsTh(kfreq)*100, 100*amorts(kfreq)]);
    end
    disp('...');
    kfreq = length(freqs);
    fprintf('mode %d: f=%.3fHz (f_th=%.3fHz, %.2f%% error), z=%.2f%%\n',...
        [kfreq, freqs(kfreq), freqsTh(kfreq), (freqs(kfreq)-freqsTh(kfreq))/freqsTh(kfreq)*100, 100*amorts(kfreq)]);
end

% deformees modales
if plotShapes
    % poles
    [defMod, poles] = eig([zeros(size(Mr)), eye(size(Mr)); -Mr\Kr, -Mr\Cr]);
    defMod = defMod(1:end/2, :);
    poles = diag(poles);
    defMod = defMod(:, imag(poles) ~= 0);
    poles = poles(imag(poles) ~= 0);
    [~, Ipoles] = sort(imag(poles));
    defMod = defMod(:, Ipoles(end/2+1:end));
    poles = poles(Ipoles(end/2+1:end));
    amorts = -real(poles)./abs(poles);
    freqs0 = abs(poles)/(2*pi);
    
    [fctDefModale, fctDefModaleAnimation] = defModalePontLin(L, linspace(0, L, N));
    
    for kf = 1:3
        f0 = freqs0(kf);
        z = amorts(kf);
        phik = defMod(:, kf);
        phik = [0; phik(1)/2; phik; phik(end)/2; 0]; % ddl manquants
        phik = phik / sqrt(phik.'*phik);
        Ik = nonPropIndex(phik);
        
        figName = sprintf('mode %u (f = %.2fHz, z = %.2f%%, I = %.2f%%)',...
            [kf, f0, 100*z, 100*Ik]);
        
%         fctDefModale(real(phik), figName);
        
        fctDefModaleAnimation(phik, figName);
        
        complexShapePlot1(phik, figName);
    end
end

end