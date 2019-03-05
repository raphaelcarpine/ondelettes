
function WvltOut=WvltComp(X,Y,WvltFreq,Qin,varargin)
p = inputParser ;
%% parametres par defaut
ZeroPaddingDef = 1;
CenterSignalDef = false;
ctDef = 5;
%%
addRequired(p,'X')
addRequired(p,'Y')
addRequired(p,'WvltFreq')
addRequired(p,'Q')
addParameter(p,'ZeroPadding',ZeroPaddingDef);
addParameter(p,'CenterSignal',CenterSignalDef);
addParameter(p,'ct',ctDef);

parse(p,X,WvltFreq,Y,Qin,varargin{:});

%
ZeroPadding = p.Results.ZeroPadding;
CenterSignal = p.Results.CenterSignal;
ct = p.Results.ct;
%%
Fs = 1/mean(X(2:end)-X(1:end-1)); %Frequence d'echantillonage

if ~iscolumn(Y)
    Y =transpose(Y); %Y en colonne
end
%%
if CenterSignal
    Y = y-mean(y);
end
%%
WvltOut = zeros(length(WvltFreq),length(Y),length(Qin));

%% Calcul
for CQ=1:length(Qin)
    Q=Qin(CQ);
    %% Choix du nb de zeros a ajouter pour chaque freq. de calcul de la CWT
    n=length(X);
    
    if ZeroPadding==0
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
            FreqSet{CSet} = N0==NSet(CSet) ; % Liste des indices des freq. associees � chaque set (auquel est associe un Nb pts pour fft)
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
        scaledftFreq = bsxfun(@times,scales,transpose(ftFreq)); % matrice ordre 2 de   a * f  avec a �chelle et f freq d'eval de fft
        
        
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