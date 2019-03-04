clear
close all
%% Infos

% Système 1ddl non-linéaire, avec composante endommageante (basée sur le
% dép. max enregistré) et composante elastoplastique. On pilote la
% participation de chaque composante avec le paramètre Chi.
% Application de l'extraction des ridges à des non-linéarités, et
% comparaison avec une estimation analytique des résultats (dans le cas d'une seule vitesse initiale non nulle).
%  Possibilité de comparer plusieurs choix de Chi en en faisant un vecteur
%  contenant plusieurs valeurs.

% Calcul de la réponse du système codé par Thuong Anh Nguyen
% http://www.theses.fr/2017PESC1032

%% Paramètres du système

temps = (0:0.0001:3); % Vecteur temps, en secondes

AccImpo = zeros(size(temps)); % AccImpo = Fext/m ; c'est la sollicitation, nulle par défaut.

Chi = [0,.5,1];% Vecteur d'éléments compris entre 0 et 1 inclus.
% Chaque élément est un ratio de transmission de la force elastoplatique
% si 1, système 100% elastoplastique ; si 0, système 100% endommageant ; mixte sinon
% Chaque element du vecteur correspond à un système indépendant étudié :
% ex : Chi = [0,.5,1] lance l'études de trois systèmes (endommageant, mixte 50-50 et elastoplastique)

AlphaP = .7 ;% Rapport entre la rigidité tangente post-seuil et la rigidité initiale
f0 = 17.0; % fréquence initiale, en Hz
xi0 = 0.02; % taux d'amortissement initial, détermine la force visqueuse.
Xe = [.0001]; % seuil d'élasticité, en mètres
u0 = 0; % condition initiale en déplacement, en mètres
v0 = .9 ; % condition initiale en vitesse, en m/s

%% Paramètres d'extraction

f_min = 13.5 ; % freq. min. de la plage de calcul
f_max = 17.5 ; % freq. max. de la plage de calcul
nb_freq = 300 ; % discrétisation de la plage de calcul
Q = 4 ; % facteur de qualité retenu

SortieIdent = 'a' ; % 'u', 'v' ou 'a' selon que l'on souhaite étudier dép., vit. ou acc.

filename='test'; % préfixe des fichiers de sortie : <filename>_Chi<Chi>_<Grandeurs>.png
dimfig=[14,16]; % [largeur,hauteur] (cm)[16,10]

%% Qqs paramètres induits
nChi = length(Chi); % Nb de systèmes différents

nDiff=strcmp('u',SortieIdent)*0 + strcmp('v',SortieIdent)*1 + strcmp('a',SortieIdent)*2; % nDiff=0,1 ou 2 si on prend u,v ou a. Utile pour calculer les courbes théoriques
m = strcmp('u',SortieIdent)*3 + strcmp('v',SortieIdent)*0 + strcmp('a',SortieIdent)*2; % Choisit le tracé A'/A (+ ou 0 ou -) f'/f fonction du type de données (dep, vit ou acc)
m = num2str(m(m>0)); % '3', '2' ou [] selon que l'on traite dep, vit ou acc
%% Calcul et identification pour chaque Chi

for CChi = 1:nChi % CChi est l'indice du système étudié
    
    [u{CChi},v{CChi},a{CChi},fint{CChi},alpha{CChi},u_p1{CChi},fint1{CChi},Y2{CChi},fint2{CChi}]=comp_mixte_06aout2017(temps,AccImpo,Chi(CChi),f0,xi0,Xe,AlphaP,u0,v0);
    
    % u, v et a : déplacement, vitesse et accélération calculés ; pour chaque Chi
    % fint, fint1, fint 2: force de rappel totale,  composante élastoplastique
    % et composante endommageante ; sans pondération
    % u_p1 : déplacement plastique
    % Y2 : variable interne d'endommagement : dép. max.
    
    %% Tracés tempo des sorties
    
    FigDVA(CChi) = figure('Name','DVA');
    subplot(2,2,1),plot(temps,u{CChi},temps,u_p1{CChi}*(Chi(CChi)>0),temps,cummax(u{CChi}),temps,[0,cumtrapz(abs(diff(u_p1{CChi})))]*(Chi(CChi)>0))
    legend('total','plas.','max(t''<t)','plas. cum.')
    leg1.Location='best';
    xlabel('time [s]')
    ylabel('disp. [m]')
    subplot(2,2,2),plot(temps,v{CChi})
    xlabel('time [s]')
    ylabel('vit. [m/s]')
    subplot(2,2,3),plot(temps,a{CChi})
    xlabel('time [s]')
    ylabel('acc. [m/s^2]')
    subplot(2,2,4),plot(u{CChi},fint{CChi})
    xlabel('dep. [m]')
    ylabel('force [N]')
    %% Extraction du ridge
    ridge{CChi} = RidgeExtract(temps,eval([SortieIdent,'{CChi}']),Q,f_min,f_max,nb_freq,'ctLeft',4,'ctRight',4);
    %% Tracés temporels des données des ridges
    FigRidge(CChi)=figure('Name','ridge');
    
    ridge0{1}=ridge{CChi}; % On sélectionne un seul ridge car on ne souhaite pas superposer les tracés temporels de systèmes différents
    
    subplot(2,2,1)
    QtyX = 'time'; %Quantité en abscisses (parmi les champs de <ridge>)
    QtyY = 'freq'; %Quantité en ordonnées
    NameX = 'time [s]'; %Nom de l'axe des abscisses
    NameY = 'freq. [Hz]'; %Nom de l'axe des ordonnées
    ScaleX = 'linear'; %'log' ou 'linear' selon l'échelle de l'axe des abscisses
    ScaleY = 'linear'; %'log' ou 'linear' selon l'échelle de l'axe des ordonnées
    LimX = [0,3]; %Plage de l'axe des abscisses
    LimY = [min(f_min),max(f_max)]; %Plage de l'axe des abscisses
    [axfreq]=RidgeQtyPlotRheo(ridge0,QtyX,QtyY,...
        NameX,NameY,ScaleX,ScaleY,LimX,LimY); % On trace avec ce qu'on a choisi précédemment
    
    subplot(2,2,2)
    QtyY = 'val';
    NameY = 'ln A';
    ScaleX = 'linear';
    ScaleY = 'log';
    LimY = [-inf,+inf];
    [axAbsLog] = RidgeQtyPlotRheo(ridge0,QtyX,QtyY,...
        NameX,NameY,ScaleX,ScaleY,LimX,LimY);
    
    subplot(2,2,3)
    QtyY = ['inv2Q',m];
    NameY = '($\sigma_f + \varsigma_f$) / f  [1]';
    ScaleX = 'linear';
    ScaleY = 'linear';
    LimY = [0,.1];
    [axdamp] = RidgeQtyPlotRheo(ridge0,QtyX,QtyY,...
        NameX,NameY,ScaleX,ScaleY,LimX,LimY);
    
    subplot(2,2,4)
    QtyY = ['bandwidth',m];
    NameY = '$\sigma_f + \varsigma_f$  [Hz]';
    ScaleX = 'linear';
    ScaleY = 'linear';
    LimY = [-inf,+inf];
    [axDiff] = RidgeQtyPlotRheo(ridge0,QtyX,QtyY,...
        NameX,NameY,ScaleX,ScaleY,LimX,LimY);
    
    linkaxes([axDiff,axdamp,axAbsLog,axfreq],'x')
end
%% Tracés temporels des données des ridges
FigRidgeAbs=figure('Name','ridge abs');

subplot(2,2,1.5)
QtyX = 'valraw';
QtyY = 'freq';
NameX = 'ln A';
NameY = 'freq. [Hz]';
ScaleX = 'log';
ScaleY = 'linear';
LimX = [-inf,+inf];
LimY = [min(f_min),max(f_max)];
[axfreqabs]=RidgeQtyPlotRheo(ridge,QtyX,QtyY,...
    NameX,NameY,ScaleX,ScaleY,LimX,LimY);
axfreqabs.ColorOrderIndex=1;

% Evaluation analytique des résultats attendus pour comparaison
for CChi=1:nChi
    legKi{CChi} = sprintf('\\chi = %.2f iden.',Chi(CChi));
    legKi{CChi+nChi} = sprintf('\\chi = %.2f theo.',Chi(CChi));
    x=logspace(log10(Xe),log10(max([u{CChi},Xe])),200);
    omega2 = Chi(CChi)*((2*pi*f0)^2*Xe./x + AlphaP*(2*pi*f0)^2*(x-Xe)./x)...
        + (1-Chi(CChi)) * ((2*pi*f0)^2*Xe./max([u{CChi},Xe]) + AlphaP*(2*pi*f0)^2*(max([u{CChi},Xe])-Xe)./max([u{CChi},Xe]));
    x=[min(abs(ridge{CChi}.valraw{1}))/((omega2(1)).^(nDiff/2)),x];
    omega2 = [omega2(1),omega2];
    hold on
    plot(log(x.*(omega2).^(nDiff/2)),sqrt(omega2)/(2*pi),'--','LineWidth',2)
end

l=legend(legKi);
l.Position(1:2)=[.7,.5];

subplot(2,2,3)
QtyY = ['inv2Q',m];
NameY = '($\sigma_f + \varsigma_f$) / f  [1]';
ScaleX = 'log';
ScaleY = 'linear';
LimY = [0,.1];
[axdampabs] = RidgeQtyPlotRheo(ridge,QtyX,QtyY,...
    NameX,NameY,ScaleX,ScaleY,LimX,LimY);
axdampabs.ColorOrderIndex=1;

% Evaluation analytique des résultats attendus pour comparaison
for CChi=1:nChi
    x=logspace(log10(Xe),log10(max([u{CChi},Xe])),200);
    omega2 = Chi(CChi)*((2*pi*f0)^2*Xe./x + AlphaP*(2*pi*f0)^2*(x-Xe)./x)...
        + (1-Chi(CChi)) * ((2*pi*f0)^2*Xe./max([u{CChi},Xe]) + AlphaP*(2*pi*f0)^2*(max([u{CChi},Xe])-Xe)./max([u{CChi},Xe]));
    sigf_f = Chi(CChi)*4*((2*pi*f0)^2-AlphaP*(2*pi*f0)^2)*Xe*(x-Xe)./(4*pi*x.^2.*omega2/2)+xi0*f0./(sqrt(omega2)/(2*pi));
    x=[min(abs(ridge{CChi}.valraw{1}))/((omega2(1)).^(nDiff/2)),x];
    omega2 = [omega2(1),omega2];
    sigf_f = [sigf_f(1),sigf_f];
    hold on
    plot(log(x.*(omega2).^(nDiff/2)),sigf_f,'--','LineWidth',2)
end


subplot(2,2,4)
QtyY = ['bandwidth',m];
NameY = '$\sigma_f + \varsigma_f$  [Hz]';
ScaleX = 'log';
ScaleY = 'linear';
LimY = [-inf,+inf];
FigName = 'Bandwidth';
[axDiffabs] = RidgeQtyPlotRheo(ridge,QtyX,QtyY,...
    NameX,NameY,ScaleX,ScaleY,LimX,LimY);
axDiffabs.ColorOrderIndex=1;

% Evaluation analytique des résultats attendus pour comparaison
for CChi=1:nChi
    x=logspace(log10(Xe),log10(max([u{CChi},Xe])),200);
    omega2 = Chi(CChi)*((2*pi*f0)^2*Xe./x + AlphaP*(2*pi*f0)^2*(x-Xe)./x)...
        + (1-Chi(CChi)) * ((2*pi*f0)^2*Xe./max([u{CChi},Xe]) + AlphaP*(2*pi*f0)^2*(max([u{CChi},Xe])-Xe)./max([u{CChi},Xe]));
    gca=axDiffabs;
    sigf_f = Chi(CChi)*8*((2*pi*f0)^2-AlphaP*(2*pi*f0)^2)*Xe*(x-Xe)./(4*pi*x.^2.*omega2)+xi0*f0./(sqrt(omega2)/(2*pi));
    sigf = sigf_f.*(sqrt(omega2)/(2*pi));
    x=[min(abs(ridge{CChi}.valraw{1}))/((omega2(1)).^(nDiff/2)),x];
    omega2 = [omega2(1),omega2];
    sigf = [sigf(1),sigf];
    hold on
    plot(log(x.*(omega2).^(nDiff/2)),sigf,'--','LineWidth',2)
end
%% Creation (si besoin) du dossier /pics/ à côté de ce .m
p=mfilename('fullpath');
p=p(1:end-length(mfilename));
p=[p,'pics'];
[~,~,~] = mkdir(p);
%% Sauvegarde des images dans le dossier /pics/
PlotList = {'FigRidge',...
    'FigRidgeAbs','FigDVA'};
for Csave = 1:length(PlotList)
    if exist(PlotList{Csave},'var')
        FigC = eval(PlotList{Csave});
        nFig = length(FigC);
        set(FigC,'Units','centimeters','Position',[0,0,dimfig])
        for Cfig=1:nFig
            filenametot = [filename,'_',sprintf('Chi%d',round(100*Chi(Cfig+(0:nChi-nFig)))),'_',PlotList{Csave}(4:end),'.png'];
            saveas(FigC(Cfig),fullfile(p,filenametot))
        end
    end
end