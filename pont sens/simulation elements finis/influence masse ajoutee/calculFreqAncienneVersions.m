clear all
close all

folder_dir = 'pont sens/simulation elements finis/resultats/testsMasseAjoutee';

%% enregistrement

enregistrement = false;

Array_Fn = [];
Array_Fn_t = [];
Array_mu = [];
Array_Lw = [];

for kSimul = 1:18
    % recherche du fichier
    name0 = sprintf('simul%d', kSimul);
    listing = dir(folder_dir);
    listingNames = {listing.name};
    fileName = [];
    for kname = 1:length(listingNames)
        name1 = strsplit(listingNames{kname}, '_');
        name1 = name1{1};
        if strcmp(name1, name0)
            fileName = listingNames{kname};
            break
        end
    end
    
    if isempty(fileName)
        error('file not found');
    end
    
    load([folder_dir, '/', fileName]);
    
    %% affichage mu, L
    
    fprintf('\n~~~~~~~~ simul %u ~~~~~~~~\n\n', kSimul);
    fprintf('mu_train = %.2f*mu_pont ; L_essieux = %.1f\n', [mu_t/mu, L_wagons]);
    
    %% temps
    
    t1 = 0; % debut d'entree du train
    t2 = ((N_bogies-1)*L_wagons + L)/c; % fin de sortie du train
    if L_wagons >= L
        t1bis = t1; % fin d'entree du train
        t2bis = t2; % debut de sortie du train
    else
        t1bis = t1 + L/c;
        t2bis = t2 - L/c;
    end
    
    fprintf('t1 = %.2f, t1'' = %.2f, t2 = %.2f, t2'' = %.2f\n', [t1, t1bis, t2, t2bis]);
    
    %% calcul freq propre, freq excitation
    try
        for kfreq = 1:1
            fprintf('freq. propre %d : %.2fHz, th. %.2fHz (%.1f%% error)\n',...
                [kfreq, freqs(kfreq), freqsTh(kfreq), (freqs(kfreq)-freqsTh(kfreq))/freqsTh(kfreq)*100]);
        end
    catch
    end
    
    fprintf('freq. excitation : %.2f\n', c/L_wagons);
    
    %% frequences pour CWT
    % phase 1 (train sur le pont)
    F_excitation = c/L_wagons;
    Fn_t = freqs(1) * sqrt(mu/(mu+mu_t)); % premier mode
    Fn2_t = freqs(2) * sqrt(mu/(mu+mu_t)); % deuxième mode
    % phase 2 (train parti du pont)
    Fn = freqs(1); % premier mode
    Fn2 = freqs(2); % deuxième mode
    
    %% resultats
    
    pos_capteurs = L/2;
    
    % wavelet capteurs
    Ycapt = getYcapt(Ytot, pos_capteurs, dx);
    
    fig = figure;
    ax = axes(fig);
    plt = plot(ax, t, Ycapt);
    
    
    %% première extraction (avec train)
    
    disp('\nfréquence avec train');
    
    lambda_t = 0.01*2*pi*Fn_t;
    
    Df = min([abs(Fn_t-F_excitation), Fn2_t-Fn_t, Fn_t]);
    Dt = 1/lambda_t;
    T = t2bis-t1bis;
    ct = 3;
    cf = 5;
    
    [Qmin, Qmax, Qz] = getBoundsQ(Fn_t, Df, Dt, T, ct, cf);
    Q = min((Qmin+Qmax)/2, Qz);
    fprintf('Qmin = %.1f, Qmax = %.1f, Qz = %.1f\n', [Qmin, Qmax, Qz]);
    fprintf('Q = %.1f\n', Q);
    
    if Qmin > Qmax || Qmin > Qz
        warning('Qmin > Qmax ou Qmin > Qz');
        continue
    end
    
    fmin = 4.3;
    fmax = 6.5;
    MaxRidges = 1;
    XLim = [t1bis, t2bis];
    fig_wvlt = WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
        'MaxRidges', MaxRidges, 'XLim', XLim);
    
    % frequence
    optionsReg = optimoptions(@lsqnonlin, 'Display', 'off');
    Fn_t = RegressionMenu('Equation', 'f', 'Param', 'f', 'Param0', Fn_t, 'OptionsRegression', optionsReg);
    
    % amort.
    Alambda_t = RegressionMenu('Equation', 'real(A)*exp(-real(l)*x)', 'Param', 'A l', 'Fit', 'log(y)',...
        'Param0', [1e-3, lambda_t], 'OptionsRegression', optionsReg);
    lambda_t = Alambda_t(2);
    
    delete(fig_wvlt);
    
    %% deuxième extraction (avec train)
    
    Df = min([abs(Fn_t-F_excitation), Fn2_t-Fn_t, Fn_t]);
    Dt = 1/lambda_t;
    T = t2bis-t1bis;
    ct = 3;
    cf = 5;
    
    [Qmin, Qmax, Qz] = getBoundsQ(Fn_t, Df, Dt, T, ct, cf);
    Q = min((Qmin+Qmax)/2, Qz);
    fprintf('Qmin = %.1f, Qmax = %.1f, Qz = %.1f\n', [Qmin, Qmax, Qz]);
    fprintf('Q = %.1f\n', Q);
    
    fmin = 4.3;
    fmax = 6.5;
    MaxRidges = 1;
    XLim = [t1bis, t2bis];
    fig_wvlt = WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
        'MaxRidges', MaxRidges, 'XLim', XLim);
    
    % frequence
    Fn_t = RegressionMenu('Equation', 'f', 'Param', 'f', 'Param0', Fn_t, 'OptionsRegression', optionsReg);
    
    % amort.
    Alambda_t = RegressionMenu('Equation', 'real(A)*exp(-real(l)*x)', 'Param', 'A l', 'Fit', 'log(y)',...
        'Param0', [1e-3, lambda_t], 'OptionsRegression', optionsReg);
    lambda_t = Alambda_t(2);
    
    delete(fig_wvlt);
    
    
    %% première extraction (sans train)
    
    disp('fréquence à vide');
    
    lambda = lambda_t;
    
    Df = min([Fn2-Fn, Fn]);
    Dt = 1/lambda;
    T = tf-t2;
    ct = 3;
    cf = 5;
    
    [Qmin, Qmax, Qz] = getBoundsQ(Fn, Df, Dt, T, ct, cf);
    Q = min((Qmin+Qmax)/2, Qz);
    fprintf('Qmin = %.1f, Qmax = %.1f, Qz = %.1f\n', [Qmin, Qmax, Qz]);
    fprintf('Q = %.1f\n', Q);
    
    fmin = 4.3;
    fmax = 6.5;
    MaxRidges = 1;
    XLim = [t2, tf];
    fig_wvlt = WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
        'MaxRidges', MaxRidges, 'XLim', XLim);
    
    % frequence
    Fn = RegressionMenu('Equation', 'f', 'Param', 'f', 'Param0', Fn, 'OptionsRegression', optionsReg);
    
    % amort.
    Alambda = RegressionMenu('Equation', 'real(A)*exp(-real(l)*x)', 'Param', 'A l', 'Fit', 'log(y)',...
        'Param0', [1e-3, lambda], 'OptionsRegression', optionsReg);
    lambda = Alambda(2);
    
    delete(fig_wvlt);
    
    %% deuxième extraction (sans train)
    
    Df = min([Fn2-Fn, Fn]);
    Dt = 1/lambda;
    T = tf-t2;
    ct = 3;
    cf = 5;
    
    [Qmin, Qmax, Qz] = getBoundsQ(Fn, Df, Dt, T, ct, cf);
    Q = min((Qmin+Qmax)/2, Qz);
    fprintf('Qmin = %.1f, Qmax = %.1f, Qz = %.1f\n', [Qmin, Qmax, Qz]);
    fprintf('Q = %.1f\n', Q);
    
    fmin = 4.3;
    fmax = 6.5;
    MaxRidges = 1;
    XLim = [t2, tf];
    WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
        'MaxRidges', MaxRidges, 'XLim', XLim);
    
    % frequence
    Fn = RegressionMenu('Equation', 'f', 'Param', 'f', 'Param0', Fn, 'OptionsRegression', optionsReg);
    
    % amort.
    Alambda = RegressionMenu('Equation', 'real(A)*exp(-real(l)*x)', 'Param', 'A l', 'Fit', 'log(y)',...
        'Param0', [1e-3, lambda], 'OptionsRegression', optionsReg);
    lambda = Alambda(2);
    
    delete(fig);
    
    close all
    
    %% calcul freqs non amorties
    
    Fn_t = abs(Fn_t + 1i*lambda_t/(2*pi));
    Fn = abs(Fn + 1i*lambda/(2*pi));
    
    fprintf('\nFn = %.2fHz (avec train), Fn = %.2fHz (sans)\n', [Fn_t, Fn]);
    
    % estimation rapport freqs
    rapport_th = sqrt(mu/(mu+mu_t));
    rapport_exp = Fn_t/Fn;
    
    fprintf('rapport de frequences :\nestimé : %.4f ; expérimental : %.4f (%.2f%% error)\n',...
        [rapport_th, rapport_exp, 100*(rapport_exp/rapport_th-1)]);
    
    %% enregistrement
    
    Array_Fn = [Array_Fn, Fn];
    Array_Fn_t = [Array_Fn_t, Fn_t];
    Array_mu = [Array_mu, mu_t/mu];
    Array_Lw = [Array_Lw, L_wagons];
    
    
end

%% enregistrement

if enregistrement
    save('pont sens/simulation elements finis/influence masse ajoutee/results', 'Array_Fn', 'Array_Fn_t', 'Array_mu', 'Array_Lw');
end