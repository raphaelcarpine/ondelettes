
function WvltOut = WvltComp(X,Y,WvltFreq,Qin,varargin)

p = inputParser ;
%% parametres par defaut
ZeroPaddingDef = true;
CenterSignalDef = false;
ctDef = 3;
%%
addRequired(p,'X')
addRequired(p,'Y')
addRequired(p,'WvltFreq')
addRequired(p,'Qin')
addParameter(p,'ZeroPadding',ZeroPaddingDef);
addParameter(p,'CenterSignal',CenterSignalDef);
addParameter(p,'ct',ctDef);

parse(p,X,WvltFreq,Y,Qin,varargin{:});

%
ZeroPadding = p.Results.ZeroPadding;
CenterSignal = p.Results.CenterSignal;
ct = p.Results.ct;
%%
if length(WvltFreq) > 500
    str = input('length(WvltFreq) > 500, continue ? ', 's');
    if ~ismember(str, {'', 'y', 'yes', 'oui', 'true'})
        error('length(WvltFreq) > 500');
    end
end

%%
Fs = 1/mean(diff(X)); %Frequence d'echantillonage

Diff = diff(X);
Diff = Diff/Diff(1);
if max(abs(Diff-1)) > 1e-5
    warning(['pas de temps non constant, erreur : ', num2str(max(abs(Diff-1)))]);
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
if max(WvltFreq) > Fs/2
    warning('Shanon');
end

%%
if CenterSignal
    Y = Y-mean(Y);
end
%%
WvltOut = zeros(length(WvltFreq),length(Y),length(Qin));

%% Calcul
for CQ=1:length(Qin)
    Q=Qin(CQ);
    %% Choix du nb de zeros a ajouter pour chaque freq. de calcul de la CWT
    n=length(X);
    
    if ~ZeroPadding
        N0=2*ceil(n/2); %multiple de 2 tout juste superieur (plus facile de gerer la periodicite de la fft)
        NbSet = 1; %un seul set de calcul
        NSet = N0; %nb de points pour fft
        FreqSet{1} = true(size(WvltFreq)); %liste des indices des freq. du set associe a N0 pts pour la fft
    else
        DeltaT = Q./(2*pi*WvltFreq); % = a * Delta t_psi
        DeltaTMaxInd = round(ct*Fs*DeltaT); % Conversion en nb de pts de ct * a * Delta t_psi = zero padding mini
        N0 = 2.^nextpow2(n+DeltaTMaxInd); % Nb pts pour fft = Puissance de 2 >= nb pts du signal + nb pts du zero-padding minimum
        NSet = unique(N0); % Liste des Nb pts pour fft uniques
        NbSet = length(NSet); % Nb de set de calcul (autant que de N0 uniques)
        FreqSet = cell(size(NSet));
        for CSet = 1:NbSet
            FreqSet{CSet} = N0==NSet(CSet) ; % Liste des indices des freq. associees à chaque set (auquel est associe un Nb pts pour fft)
        end
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
        
        ftFreq=(0:N/2-1)*Fs/N; %Freq d'eval de la fft, f positif
        
        n_cau = (2*Q.^2 - 1/2); % parametre de l'ondelette de Cauchy
        fCau = (n_cau)/(2*pi);% freq. centrale de l'ondelette
        scales = fCau./WvltFreq1 ; % echelles associees aux freq. de calcul
        scaledftFreq = bsxfun(@times,scales,transpose(ftFreq)); % matrice ordre 2 de   a * f  avec a échelle et f freq d'eval de fft
        
        
        ftWvlt{CSet} = exp(n_cau-2*pi*scaledftFreq + n_cau*log(scaledftFreq/fCau)); % Calcul Ondelettes cauchy en a * f pour f>0
        
        WvltNorm = 1/2;    % Coef de normalisation pour avoir max(ftWvlt) = 2
        
        ftWvlt{CSet} = cat(1,ftWvlt{CSet},zeros(N/2,length(WvltFreq1))); % on rajoute des zeros pour a*f<0
    end
    %% integrande
    % produit de la fft des ondelettes scalees avec la fft du signal
    integrande = cell(size(NSet));
    for CSet = 1:NbSet
        integrande{CSet}=bsxfun(@times,ftWvlt{CSet},ftSignal{CSet}); % Calcul par set
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
    WvltOut(:,:,CQ)=transpose(wavelet1)/WvltNorm; % On termine avec le coef de normalisation
end