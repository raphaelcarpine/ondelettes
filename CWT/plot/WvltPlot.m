function WvltPlot(X,Y,WvltFreq,Q,varargin)%freq_scale
p = inputParser ;

%% Param def
WvltNameDef = 'Cauchy'; % famille d'ondelettes, seul Cauchy ici 
ZeroPaddingDef = 1; % zero-padding par defaut
ctDef = 3; % ct=5 par defaut (pour zero-padding ET effets de bord)
PlotScaleDef = 'log'; % par defaut, on trace la CWT en module
FreqScaleDef = 'lin'; % echelle de l'axe des frequences lineaire par defaut
TitleDef = []; % titre de la figure
VisibleDef = 'on'; % figure visible par defaut (si bcp de pts, figure lourde a afficher... 'off' utile pour sauvegarder l'image sans l'afficher)
%%
addRequired(p,'X')
addRequired(p,'Y')
addRequired(p,'WvltFreq')
addRequired(p,'Q')
addParameter(p,'ZeroPadding',ZeroPaddingDef);
addParameter(p,'ctEdgeEffects',ctDef);
addParameter(p,'PlotScale',PlotScaleDef);

addParameter(p,'FreqScale',FreqScaleDef);
addParameter(p,'Title',TitleDef);
addParameter(p,'Visible',VisibleDef);


parse(p,X,Y,WvltFreq,Q,varargin{:});


ctEdgeEffects = p.Results.ctEdgeEffects;

WvltName = WvltNameDef;
ZeroPadding = p.Results.ZeroPadding;
PlotScale = p.Results.PlotScale;
FreqScale = p.Results.FreqScale;
Title = p.Results.Title;
Visible = p.Results.Visible;


%% Creation de la figure (avec ou sans titre)
if isempty(Title)
    figure('Visible',Visible)
else
    figure('Name',Title,'NumberTitle','off','Visible',Visible)
end
%% Calcul de la CWT
[wavelet, ctZeroPadding] = WvltComp(X,Y,WvltFreq,Q, 'ZeroPadding',ZeroPadding);
%% on veut Y colonne
if ~iscolumn(Y)
    Y=transpose(Y);
end
%% Freq. d'echantillonage
Fs=1/mean(X(2:end)-X(1:end-1));
%% plot temporel
subplot(2,2,1)
plot(X,Y)
ax1=gca;
xlabel('time [T]')
ylabel('signal [S]')
%% plot CWT fonction de l'echelle choisie
subplot(2,2,3)
if strcmp(PlotScale,'log')
    pcolor(X,WvltFreq,log10(abs(wavelet)));
    ColorbarLabel='log_{10} |CWT|';
elseif strcmp(PlotScale,'abs')
    pcolor(X,WvltFreq,abs(wavelet));
    ColorbarLabel= '|CWT|';
elseif strcmp(PlotScale,'squ')
    pcolor(X,WvltFreq,abs(wavelet).^2);
    ColorbarLabel='|CWT|^{-1/2}';
end
c=colorbar('eastoutside');
c.Title.String = ColorbarLabel;
c.Position=[.8,.6,.05,.25];

ax2=gca;
xlabel('time [T]')
ylabel('frequency [T]^{-1}')
shading flat
colormap(jet)
hold on
%% zones d'effets de bord
if ctEdgeEffects>0
    delta_t = ctEdgeEffects*Q./(2*pi*WvltFreq); % = a * Delta t_psi, en considerant que 2*pi*Delta f_psi * Delta t_psi = 1/2
    
    range1 = delta_t<((X(end)-X(1))/2); % liste des indices des freq. telles que les effets de bord droite et gauche ne se superposent pas 
    delta_t1 = delta_t(range1); % liste des delta_t tels que effets de bord droite et gauche ne se superposent pas
    
    % Creation d'un "patch" pour colorer la zone d'effets de bord. On
    % cherche simplement les points qui delimite la surface (On a fait
    % attention a traiter correctement les recouvrements d'effets de
    % bord... Mais c'est un detail de tracage)
    Xedge = [X(1),X(1),X(1)+flip(delta_t1)];
    Yedge = [WvltFreq(1),WvltFreq(end),flip(WvltFreq(range1))];
    Xedge = [Xedge,(X(1)+X(end))/2];
    Yedge = [Yedge,max(WvltFreq(1),ctEdgeEffects*Q/(pi*(X(end)-X(1))))];
    Xedge = [Xedge,  X(end)-delta_t1,X(end),X(end)];
    Yedge = [Yedge, WvltFreq(range1),WvltFreq(end),WvltFreq(1)];
    
    patch(Xedge,Yedge,'white','EdgeColor','none','FaceAlpha',.35);
    
    plot(X(1)+delta_t,WvltFreq,'black','LineWidth',1.5,'LineStyle',':') % limite effets de bord gauche
    plot(X(end)-delta_t,WvltFreq,'black','LineWidth',1.5,'LineStyle',':') % limite effets de bord droite
end
%% Plot de la fft
df = min(WvltFreq(2:end)-WvltFreq(1:end-1)); % df = pas minimum des frequences de calcul de la CWT
N = 2^nextpow2(max([Fs/df,length(X)])); % On choisit N = nb de pts pour fft tel que le pas freq. soit au moins inf. a df
freq=(0:N/2)*Fs/N; % Liste des freq. (positives)

ft_Y_w0=fft(Y,N,1)/Fs; % Transformee de Fourier

esdY_w = abs(ft_Y_w0(1:N/2+1)).^2; % module de la TF au carré

range_f = freq>=WvltFreq(1) & freq<=WvltFreq(end); % On limite le plot aux limites des freq. de calcul.

freq=freq(range_f); % On retire donc les freq. hors freq. de calcul de CWT
esdY_w=esdY_w(range_f); % et on retire donc l'ESD hors freq. de calcul de CWT

subplot(2,2,4)
hold on
plot(esdY_w,freq,'DisplayName','ESD')
%% Estimation de l'ESD par la CWT
if ZeroPadding == 1 % La relation est valide sur CWT non coupee : si zeropadding, alors la CWT a ete coupee apres calcul... il faut la (re)calculer non coupee (sans zeropadding)
    wavelet=WvltComp(X,Y,WvltFreq,Q, 'ZeroPadding',0);
end
nCau = (2*Q.^2 - 1/2); % parametre de l'ondelette
Cpsi = exp((1-nCau)*log(4) + 2*nCau -2*nCau*log(nCau) + gammaln(2*nCau)); % Constante d'admissibilite de l'ondelette de param nCau
WvltESD = (1/Cpsi)*sum(abs(wavelet).^2,2)./(Fs*transpose(WvltFreq)); % Calcul de ESD_psi = int(|CWT|^2 db) /(Cpsi * freq) et 2*int(ESD_psi,0,+inf) = int(Yin^2,-inf,+inf)  - - - comparable à 2*int(ESD,0,+inf) = int(Yin^2,-inf,+inf)
plot(WvltESD,WvltFreq,'LineWidth',1.5,'DisplayName','$\mathrm{ESD}_\psi$') % plot de ESD_psi

%% Axes ESD
ax3=gca;

xlabel('[S]^2[T]^2')
ylabel('frequency [T]^{-1}')

lgd = legend('show');
lgd.Location = 'best';
lgd.Interpreter = 'latex';
%% linkaxes
ax1.XLim = [X(1),X(end)];
ax3.YLim = [WvltFreq(1),WvltFreq(end)];
linkaxes([ax2,ax1],'x');
linkaxes([ax2,ax3],'y');
%% echelle axe freq.
ax2.YScale = FreqScale;
ax3.YScale = FreqScale;
    
%% Titre et parametres
subplot(2, 2, 2);
title(Title);
axis off
ParamX = -.1;

text(ParamX,.5,sprintf('Q_\\psi = %3.1f',Q));

text(ParamX,.8,sprintf('\\psi: %s wavelet',WvltName));
if ZeroPadding==1
    text(ParamX,.2,sprintf('\\itc_{t,\\rmzp}\\rm \\geq %.0f', ctZeroPadding));
else
    text(ParamX,.2,sprintf('\\itc_{t,\\rmzp}\\rm = %.0f',0));
end
text(ParamX,-.1,sprintf('\\itc_{t,\\rmef}\\rm = %.0f', ctEdgeEffects));