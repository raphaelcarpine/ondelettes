function ridge = RidgeExtract(X,Y,Q,fmin,fmax,NbFreq,varargin)
%%
p = inputParser ;
%% parametres par defaut
etDef = 1;
efDef = 1;
ctLeftDef = 3;
ctRightDef = 3;
% ctLeftDef = 0;
% ctRightDef = 0;
NbMaxRidgesDef = 100;
NbMaxParallelRidgesDef = 1;
MinModuDef = 0;
LengthMinRidgeDef = 0;
ReleaseTimeDef = X(1);

%%
addRequired(p,'X')
addRequired(p,'Y')
addRequired(p,'Q')
addRequired(p,'fmin')
addRequired(p,'fmax')
addRequired(p,'NbFreq')
addParameter(p,'et',etDef);
addParameter(p,'ef',efDef);
addParameter(p,'ctLeft',ctLeftDef);
addParameter(p,'ctRight',ctRightDef);
addParameter(p,'NbMaxRidges',NbMaxRidgesDef);
addParameter(p,'NbMaxParallelRidges',NbMaxParallelRidgesDef);
addParameter(p,'MinModu',MinModuDef);
addParameter(p,'LengthMinRidge',LengthMinRidgeDef);
addParameter(p,'ReleaseTime',ReleaseTimeDef);

parse(p,X,Y,Q,fmin,fmax,NbFreq,varargin{:});

%
et = p.Results.et ;
ef = p.Results.ef ;
ctLeft =  p.Results.ctLeft;
ctRight =  p.Results.ctRight;
NbMaxRidges = p.Results.NbMaxRidges;
NbMaxParallelRidges = p.Results.NbMaxParallelRidges;
MinModu = p.Results.MinModu;
LengthMinRidge = p.Results.LengthMinRidge;
ReleaseTime = p.Results.ReleaseTime;
%%
if iscolumn(X)
    X=transpose(X);
end
%%  %%  %%  %%  %%  %%

WvltFreq = linspace(fmin,fmax,NbFreq); % Freq de calcul de la CWT

wavelet= WvltComp(X,Y,WvltFreq,Q); % Calcul CWT
% wavelet = WvltCompTemp(X, Y, WvltFreq, length(X), Q);

mesu = raindrop(wavelet); % Recherche des maximum locaux en echelle
mesu = abs(mesu);

mesu(mesu<MinModu) = NaN; % retrait des max. locaux < MinModu

%retrait des max locaux paralleles excessifs
for iT = 1:length(X)
    localMax = mesu(:, iT);
    localMax = localMax(~isnan(localMax));
    if length(localMax) > NbMaxParallelRidges
        localMax = sort(localMax, 'descend');
        localMin = localMax(NbMaxParallelRidges);
        for ifreq = 1:NbFreq
            if mesu(ifreq, iT) < localMin
                mesu(ifreq, iT) = NaN;
            end
        end
    end
end

CauN = (2*Q^2 - 1/2); % param ondelette
CauFreq = (CauN)/(2*pi); %freq centrale ondelette
DeltaTPsi = 1/sqrt(2*CauN-1); % Delta T_psi

Scales = CauFreq./WvltFreq; % echelles a
DeltaTimeLeft = ctLeft*Scales*DeltaTPsi; % ct * a * Delta T_psi gauche = largeur effet de bord gauche
DeltaTimeRight = ctRight*Scales*DeltaTPsi; % ct * a * Delta T_psi droite = largeur effet de bord droite
RangeLeft = bsxfun(@minus,X,transpose(DeltaTimeLeft))>X(1); % Zone d'effets de bord gauche
RangeRight = bsxfun(@plus,X,transpose(DeltaTimeRight))<X(end); % Zone d'effets de bord droite
mesuEdge = mesu;
mesuEdge(~(RangeLeft&RangeRight))=NaN; % Max locaux hors zone effets de bord seulement
%% si Y_in = 0 au debut ou fin : zero-padding non voulu. On le retire
IndBegin = find(Y~=0,1); % eventuel zero-padding ou apparente (succession de 0 au début du signal)
IndEnd = find(Y~=0,1,'last');% Idem fin de signal
%%
Fs=1/mean(X(2:end)-X(1:end-1)); %Freq echantillonnage
LengthMin = max(3,round(LengthMinRidge*Fs)); % Nb points mini d'un ridge
if iscolumn(X)
    X=transpose(X); % X_in en ligne
end
%%
[ny,nx]=size(mesu); % ny = nb de freq de CWT, nx = nb de points du signal

%%
ridge.time = {};
ridge.val = {};
ridge.freq = {};
for C_r=1:NbMaxRidges % Pour chaque ridge
    [M2,I] = max(mesuEdge(:)); % Maximum global hors effets de bord pour initialiser le chainage
    
    if M2<MinModu || isnan(M2)
        break
    end
    
    ind_freq = NaN(size(X)) ; % Indice des fréquences du ridge
    ridge.time{C_r} = NaN(size(X)) ; % Instants associés
    ridge.val{C_r}  = NaN(size(X)) ; % Valeur (complexe) de la CWT
    
    [a,b] = ind2sub(size(mesu),I);    % Conversion de l'index du max global en freq,instant
    
    ind_freq(1) = a ;
    ridge.time{C_r}(1) = b ;
    ridge.val{C_r}(1)  = wavelet(a,b) ;
    
    C_ind = 2;
    while (ridge.time{C_r}(C_ind-1)+et<=nx) % On chaine vers l'avant
        
        [M1,I1] = max(mesu(bound(ind_freq(C_ind-1)+(-ef:ef),1,ny),bound(ridge.time{C_r}(C_ind-1)+(1:et),1,nx))); % Max par colonne
        [M2,I2] = max(M1); % Max des max par colonne
        
        [a,b] = ind2sub([ef*2+1,et],I1(I2)) ; % Conversion en indices
        
        I_f = bound(a + ind_freq(C_ind-1) - ef - 1,1,ny) ; %Correction des indices de freq.
        I_t = b + ridge.time{C_r}(C_ind-1) ; % Correction des indices d'instant
        
        if (M2 <= MinModu)||(isnan(M2))||(I_f==1)||(I_f==ny) % Condition de fin : module trop petit ou extremite de la fenêtre frequentielle choisie
            break
        else
            ind_freq(C_ind) = I_f ; % Assignation freq.
            ridge.time{C_r}(C_ind) = I_t ; % Assignation instant
            ridge.val{C_r}(C_ind)  = wavelet(I_f,I_t) ; % Assignation valeur
            
            C_ind = C_ind+1; % Increment
        end
    end
    NoNaN = ~isnan(ridge.time{C_r}); % On liste les instants ou on a identifier le ridge jusqu'ici
    ind_freq(NoNaN) = flip(ind_freq(NoNaN) ); % On en inverse le sens (pour chainer vers l'arriere)
    ridge.time{C_r}(NoNaN) = flip(ridge.time{C_r}(NoNaN)); % idem instants
    ridge.val{C_r}(NoNaN) = flip(ridge.val{C_r}(NoNaN)); % idem valeurs
    
    while ridge.time{C_r}(C_ind-1)-et>=1 %On chaine vers l'arrière
        
        [M1,I1] = max(mesu(bound(ind_freq(C_ind-1)+(-ef:ef),1,ny),bound(ridge.time{C_r}(C_ind-1)-(1:et),1,nx))); % Max par colonne
        [M2,I2] = max(M1); % Max des max par colonne
        
        [a,b] = ind2sub([ef*2+1,et],I1(I2)) ; % Conversion en indices
        
        I_f = bound(a + ind_freq(C_ind-1) - ef - 1,1,ny) ; %Correction des indices de freq.
        I_t = -b + ridge.time{C_r}(C_ind-1) ; % Correction des indices d'instant
        
        if (M2 <= MinModu)||(isnan(M2))||(I_f==1)||(I_f==ny) % Condition de fin : module trop petit ou extremite de la fenêtre frequentielle choisie
            break
        else
            ind_freq(C_ind) = I_f ; % Assignation freq.
            ridge.time{C_r}(C_ind) = I_t ; % Assignation instant
            ridge.val{C_r}(C_ind)  = wavelet(I_f,I_t) ; % Assignation valeur
            
            C_ind = C_ind+1; %Increment (comme on a inverser le sens, on incremente bien en +1 pour aller vers l'arriere
        end
    end
    
    NoNaN = 1:C_ind-1; % On liste les instants retenus
    ridge.time{C_r}=ridge.time{C_r}(NoNaN); % On retire les instants sans ridge
    ridge.val{C_r}=ridge.val{C_r}(NoNaN); % idem valeur
    ind_freq = ind_freq(NoNaN); % idem freq
    
    [ridge.time{C_r},Isort] = sort(ridge.time{C_r}); % On remet les instants dans l'ordre croissant
    ind_freq = ind_freq (Isort); % Idem freq (ordre chronologique)
    ridge.val{C_r} = ridge.val{C_r}(Isort); % Idem valeur (ordre chronologique)
    
    mesu(sub2ind([ny,nx],ind_freq,ridge.time{C_r})) = NaN; % On retire le ridge de nos maximum locaux potentiellement associés a un ridge pour ne pas retomber toujours sur le meme
    mesuEdge(sub2ind([ny,nx],ind_freq,ridge.time{C_r})) = NaN; % Idem hors effets de bord
    
    % On convertit les indices en frequences
    WvltFreqLocal = [WvltFreq(ind_freq-1); WvltFreq(ind_freq); WvltFreq(ind_freq+1)];
    ampl = abs(wavelet);
    mesuEdgeLocal = NaN(3, length(ind_freq));
    for iT = 1:length(ind_freq)
        mesuEdgeLocal(:,iT) = [...
            ampl(ind_freq(iT)-1, ridge.time{C_r}(iT));...
            ampl(ind_freq(iT), ridge.time{C_r}(iT));...
            ampl(ind_freq(iT)+1, ridge.time{C_r}(iT))];
    end
    ridge.freq{C_r} = localMax3Points(WvltFreqLocal, mesuEdgeLocal);
    
    % On convertit les indices en temps
    ridge.time{C_r} = X(ridge.time{C_r});
end
%% retrait des ridges trop courts
d_r=0;
for C_r=1:length(ridge.time)
    if length(ridge.time{C_r-d_r})<=LengthMin        
        FieldList = fieldnames(ridge);
        for C_Field = 1:length(FieldList)
            FieldName = FieldList{C_Field};
            ridge.(FieldName)=ridge.(FieldName)([1:C_r-d_r-1,C_r-d_r+1:length(ridge.(FieldName))]);
        end
        
        d_r=d_r+1;
    end
end

if isempty(ridge.time)
    return
end
%% Calcul de la phase (argument complexe)
for C_r = 1:length(ridge.time)
    ridge.pha{C_r} = angle(ridge.val{C_r});
end
%% frequence instantanee alternative plus lisse, en reliant les discontinuites (saut d'une freq. a une autre)
for C_r = 1:length(ridge.time)
    FreqLeft = ridge.freq{C_r}(1:end-1);
    FreqRight = ridge.freq{C_r}(2:end);
    Range = diff(ridge.freq{C_r})~=0; % Liste des sauts de freq.
    Range(end)=1;
    ListFreq = [FreqLeft(1),(  FreqLeft(Range) + FreqRight(Range)  )/2]; % a chaque saut, on associe la moyenne de la freq avant et après
    ListTime = ridge.time{C_r}(logical([1,Range])); % instants associes aux sauts de freq.
    ridge.freq2{C_r} = interp1(ListTime,ListFreq,ridge.time{C_r}); % interpolation avec ces valeurs pour tous les autres instants
    
    ridge.freq2{C_r}(1)=ridge.freq{C_r}(1); % correction premier point
    if isnan(ridge.freq2{C_r}(end))
        ridge.freq2{C_r}(end)=ridge.freq{C_r}(end); % correction dernier point
    end
end
%% evaluation de la dissipation
for C_r = 1:length(ridge.time)
    % Si mesure de vitesse
    ridge.bandwidth{C_r} = -diff(log(abs(ridge.val{C_r})))*Fs; % -A'/A
    ridge.bandwidth{C_r} = [ridge.bandwidth{C_r}(1),ridge.bandwidth{C_r}]/(2*pi); % = -A'/2*pi*A
    ridge.inv2Q{C_r} = ridge.bandwidth{C_r}./(ridge.freq2{C_r}); % = -A'/2*pi*f*A
    
    % Si mesure d'acceleration et frequence variable
    ridge.bandwidth2{C_r} = -diff(log(abs(ridge.val{C_r}./ridge.freq2{C_r})))*Fs; % -A'/A + f'/f
    ridge.bandwidth2{C_r} = [ridge.bandwidth2{C_r}(1),ridge.bandwidth2{C_r}]/(2*pi); % -A'/2*pi*A + f'/2*pi*f
    ridge.inv2Q2{C_r} = ridge.bandwidth2{C_r}./(ridge.freq2{C_r}); % -A'/2*pi*A*f + f'/2*pi*f^2
    
    % Si mesure de deplacement et frequence variable
    ridge.bandwidth3{C_r} = -diff(log(abs(ridge.val{C_r}.*ridge.freq2{C_r})))*Fs; % -A'/A - f'/f
    ridge.bandwidth3{C_r} = [ridge.bandwidth3{C_r}(1),ridge.bandwidth3{C_r}]/(2*pi); % -A'/2*pi*A - f'/2*pi*f
    ridge.inv2Q3{C_r} = ridge.bandwidth3{C_r}./(ridge.freq2{C_r}); % -A'/2*pi*A*f + f'/2*pi*f^2
    
end

%% NaN pour effets de bord, etc.
FieldList = fieldnames(ridge);
for C_Field = 1:length(FieldList)
    FieldName = FieldList{C_Field};
    if ~strcmp(FieldName,'time')
        ridge.([FieldName,'raw']) = ridge.(FieldName); % on nomme "ridge.trucraw" la grandeur "truc" hors et dans effets de bord
    end
end
FieldList = fieldnames(ridge);
for C_r = 1:length(ridge.time)
    Scales = CauFreq./ridge.freq{C_r}; % Echelles associees au ridge
    DeltaTimeLeft = ctLeft*Scales*DeltaTPsi; % ct * Delta t_psi/a associes a gauche
    DeltaTimeRight = ctRight*Scales*DeltaTPsi; % ct * Delta t_psi/a associes a droite
    
    
    RangeEdgeLeft = (ridge.time{C_r}-DeltaTimeLeft)>=X(IndBegin); % Liste des instants ou t + ct*Delta t_psi/a reste dans les limites du signal (gauche)
    RangeEdgeRight = (ridge.time{C_r}+DeltaTimeRight)<=X(IndEnd); % Liste des instants ou t + ct*Delta t_psi/a reste dans les limites du signal (droite)
    RangeEdge = RangeEdgeLeft & RangeEdgeRight; % intersection des deux listes = liste des instants hors effets de bord pour le ridge considere
    
    RangeRelease = (ridge.time{C_r}-DeltaTimeLeft)>=ReleaseTime; % Liste des instants verifiant la meme condition a partir de l'instant de relachement (debut de la reponse libre)
    
    RangeFinal = RangeEdge & RangeRelease; % intersection des listes
    
    %%
    for C_Field = 1:length(FieldList)
        FieldName = FieldList{C_Field};
        if ~strcmp(FieldName,'time') && ~strcmp(FieldName(end-2:end),'raw')
            ridge.(FieldName){C_r}(~RangeFinal) = NaN;  % NaN pour les pts proches des  bord et/ou proches du relachement
        end
    end
    
end
%% Retrait des ridges entierement dans la zone d'effets de bord
d_r=0;
for C_r=1:length(ridge.time)
    if all(isnan(ridge.val{C_r-d_r}))
        FieldList = fieldnames(ridge);
        
        for C_Field = 1:length(FieldList)
            FieldName = FieldList{C_Field};
            ridge.(FieldName) = ridge.(FieldName)([1:C_r-d_r-1,C_r-d_r+1:length(ridge.(FieldName))]);
        end
        
        d_r=d_r+1;
    end
end