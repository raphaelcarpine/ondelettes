function [WvltOut, ctZeroPadding] = WvltComp(X,Y,WvltFreq, Q, varargin)

p = inputParser ;
%% parametres par defaut
ZeroPaddingDef = true;
CenterSignalDef = false;
MotherWaveletDef = 'cauchy'; %'littlewood-paley';
ctDef = 3;
MeanOverFreqFuncDef = []; % mode mean(func(wvlt), 1)
XindexOutDef = []; % mode WvltOut = Wvlt(:, XindexOut);
MeanSquareOverXelementsDef = {};

%%
addRequired(p,'X')
addRequired(p,'Y')
addRequired(p,'WvltFreq')
addRequired(p,'Q')
addParameter(p,'ZeroPadding',ZeroPaddingDef);
addParameter(p,'CenterSignal',CenterSignalDef);
addParameter(p,'MotherWavelet',MotherWaveletDef);
addParameter(p,'ct',ctDef);
addParameter(p,'MeanOverFreqFunc', MeanOverFreqFuncDef);
addParameter(p,'XindexOut', XindexOutDef);
addParameter(p,'MeanSquareOverXelements', MeanSquareOverXelementsDef);

parse(p,X,WvltFreq,Y,Q,varargin{:});

%
ZeroPadding = p.Results.ZeroPadding;
CenterSignal = p.Results.CenterSignal;
MotherWavelet = p.Results.MotherWavelet;
ct = p.Results.ct;
MeanOverFreqFunc = p.Results.MeanOverFreqFunc;
XindexOut = p.Results.XindexOut;
MeanSquareOverXelements = p.Results.MeanSquareOverXelements;

ct = ct * ZeroPadding;

%% 
if islogical(XindexOut) || length(XindexOut) == length(X)
    XindexOut = find(XindexOut); % conversion de [0, 1, 1, 0, 1] en [2, 3, 5]
end

%%
if length(WvltFreq) > 500
    str = input('length(WvltFreq) > 500, continue? (y/n) ', 's');
    if ~ismember(str, {'', 'y', 'yes', 'o', 'oui', 'true'})
        error('length(WvltFreq) > 500');
    end
end

%% array size
Fs = 1/mean(diff(X)); %Frequence d'echantillonage

Diff = diff(X);
Diff = Diff/Diff(1);
if max(abs(Diff-1)) > 1e-5
    warning(['non-constant time step, error : ', num2str(max(abs(Diff-1)))]);
end

% vecteurs colone
if size(X, 2) ~= 1
    X = transpose(X);
end
if size(Y, 2) ~= 1
    Y = transpose(Y);
end

if size(X, 2) ~= 1 || size(Y, 2) ~= 1 || ~isequal(size(X), size(Y))
    error(['array size problem (', num2str(size(X)), ' & ', num2str(size(Y)), ')']);
end

%% shanon
if max(WvltFreq) / (Fs/2) > 1 + 1e-6
    warning('Shanon');
end

%%
if CenterSignal
    Y = Y-mean(Y);
end

%% signal long (waitbar)

longSignal = 1e-9 * length(X)*log( length(X))*  length(WvltFreq) > 1;

longSignal = longSignal && isempty(MeanOverFreqFunc) && isempty(XindexOut);

%% Calcul

[FTpsi, DeltaTfunc] = FTpsi_DeltaT(Q, MotherWavelet);

if isempty(MeanOverFreqFunc) && isempty(XindexOut) && ~longSignal % mode normal
    WvltFreqTot = {WvltFreq};
else % mode moyenne sur les fréquences ou temps particuliers (économie de la ram)
    WvltFreqTot = num2cell(WvltFreq);
    if ~isempty(MeanOverFreqFunc) && ~isempty(XindexOut)
        WvltOutTot = zeros(1, length(XindexOut));
    elseif ~isempty(MeanOverFreqFunc)
        WvltOutTot = zeros(1, length(X));
    elseif ~isempty(XindexOut)
        WvltOutTot = nan(length(WvltFreq), length(XindexOut));
    elseif longSignal
        WvltOutTot = nan(length(WvltFreq), length(X));
    end
end

% mean square over X
MeanSquareOverX = nan(0, length(MeanSquareOverXelements));

% waitbar
if ~isempty(MeanOverFreqFunc) || ~isempty(XindexOut) || longSignal
    [initWaitBar, updateWaitBar, closeWaitBar] = getWaitBar(length(WvltFreqTot),...
        'displayTime', 0, 'windowTitle', 'Computing CWT');
    initWaitBar();
end

for k_wvlt_freq = 1:length(WvltFreqTot)
    WvltFreq = WvltFreqTot{k_wvlt_freq};
    
    %% Choix du nb de zeros a ajouter pour chaque freq. de calcul de la CWT
    n = length(X);
    
    DeltaT = max(DeltaTfunc(WvltFreq), Q./(2*pi*WvltFreq)); % = a * Delta t_psi, on prend le max avec le cas mu_psi=1/2 pour l'ondelette harmonique
    DeltaTMaxInd = round(ct*Fs*DeltaT); % Conversion en nb de pts de ct * a * Delta t_psi = zero padding mini
    N0 = 2.^nextpow2(n+DeltaTMaxInd); % Nb pts pour fft = Puissance de 2 >= nb pts du signal + nb pts du zero-padding minimum
    NSet = unique(N0); % Liste des Nb pts pour fft uniques
    NbSet = length(NSet); % Nb de set de calcul (autant que de N0 uniques)
    FreqSet = cell(size(NSet));
    for CSet = 1:NbSet
        FreqSet{CSet} = N0==NSet(CSet) ; % Liste des indices des freq. associees à chaque set (auquel est associe un Nb pts pour fft)
    end
    %% fft signal
    ftSignal = cell(size(NSet));
    for CSet = 1:NbSet
        ftSignal{CSet} =fft(Y,NSet(CSet)); % Calcul fft signal avec nb de pts determine avant pour chaque set
    end
    %% wavelet
    ftWvlt = cell(size(NSet));
    for CSet = 1:NbSet
        N = NSet(CSet); % Nb pts fft associe au set
        WvltFreq1 = WvltFreq(FreqSet{CSet}); % freq. associees au set
        
        ftFreq = (0:N/2-1)*Fs/N; %Freq d'eval de la fft, f positif
        ftFreqNeg = (-N/2+1:-1)*Fs/N; %Freq d'eval de la fft, f negatif
        scales = 1./WvltFreq1 ; % echelles associees aux freq. de calcul (a = f_psi^0/a = 1/a)
        scaledftFreq = transpose(ftFreq) * scales; % matrice ordre 2 de   a * f  avec a échelle et f freq d'eval de fft
        scaledftFreqNeg = transpose(ftFreqNeg) * scales;
        
        ftWvlt{CSet} = [FTpsi(scaledftFreq);... % Calcul Ondelettes Morlet en a * f
            zeros(1, length(WvltFreq1));... % freq en N/2 inutilisable car N pair
            FTpsi(scaledftFreqNeg)]; % Calcul Ondelettes Morlet en a * f
        
        WvltNorm = 1/2;    % Coef de normalisation pour avoir max(ftWvlt) = 2
    end
    %% integrande
    % produit de la fft des ondelettes scalees avec la fft du signal
    integrande = cell(size(NSet));
    for CSet = 1:NbSet
        integrande{CSet} = ftWvlt{CSet} .* ftSignal{CSet}; % Calcul par set
    end
    %% Sortie
    % On a ajoute des zeros en augmentant le nb de points pour la fft, on
    % ajoute donc des points de calcul de la CWT
    % On decide alors de restreindre la longueur de la CWT au signal initial
    
    % On calcule la cwt
    wavelet0 = cell(size(NSet));
    for CSet = 1:NbSet
        wavelet0{CSet} = ifft(integrande{CSet},[],1); % Calcul CWT
        wavelet0{CSet} = wavelet0{CSet}(1:n,:); % Restriction signal initial
    end
    
    wavelet1=zeros(length(X),length(WvltFreq));
    for CSet = 1:NbSet
        wavelet1(:,FreqSet{CSet}) = wavelet0{CSet}; % On agrege le calcul de chaque set
    end
    WvltOut = transpose(wavelet1)/WvltNorm; % On termine avec le coef de normalisation
    
    %% moyennage sur t
    MeanSquareOverXnextLines = [];
    for k_el = 1:length(MeanSquareOverXelements)
        MeanSquareOverXnextLines = [MeanSquareOverXnextLines,...
            mean(abs(WvltOut(:, MeanSquareOverXelements{k_el})).^2, 2)];
    end
    MeanSquareOverX = [MeanSquareOverX; MeanSquareOverXnextLines];
    
    %% mode moyenne sur les fréquences ou temps particuliers (économie de la ram)
    
    if ~isempty(MeanOverFreqFunc) && ~isempty(XindexOut)
        WvltOutTot = WvltOutTot + MeanOverFreqFunc(WvltOut(:, XindexOut));
    elseif ~isempty(MeanOverFreqFunc)
        WvltOutTot = WvltOutTot + MeanOverFreqFunc(WvltOut);
    elseif ~isempty(XindexOut)
        WvltOutTot(k_wvlt_freq, :) = WvltOut(:, XindexOut);
    elseif longSignal
        WvltOutTot(k_wvlt_freq, :) = WvltOut;
    else % mode normal
        WvltOutTot = WvltOut;
    end
    
    %% waitbar
    
    if ~isempty(MeanOverFreqFunc) || ~isempty(XindexOut) || longSignal
        updateWaitBar(k_wvlt_freq);
    end
    
end

% waitbar
if ~isempty(MeanOverFreqFunc) || ~isempty(XindexOut) || longSignal
    closeWaitBar();
end

%% moyennage sur t
WvltOutTot = [WvltOutTot, MeanSquareOverX];


%%
WvltOut = WvltOutTot;

%%
ctZeroPadding = ct; % ZeroPadding*ct

end