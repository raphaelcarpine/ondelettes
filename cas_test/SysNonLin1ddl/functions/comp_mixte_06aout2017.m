function[u,v,a,fint,alpha,u_p1,fint1,Y2,fint2]=comp_mixte_06aout2017(t_cal,acc_cal,ki,f0,xi_ini,XE,alpha_p,u_0,v_0)
% le 20 fev 2016
% ce programme est le comportement mixte
% prendre en compte le comportement elasto plqstique
% et le modele de Masars pour le comportement d'endommagement
% le 31/10/2016: correction faute forces composantes

format long

% Chargement 
% t_cal doit être une colonne
% acc_cal doit être une ligne
% Caracteristiques of oscillateur 
% ki : facteur de preponderante 2 comportements
% ki = 1: elastoplastique
% ki = 0: emdommagement
w_ini=2*pi*f0; % pulsation propre de l'oscillateur
% xi_ini=xi; % initial damping of the oscillateur
% XE =0.05; % seuil d'elasticite en deplacement
alpha_e = 1; % pour calculer le facteur de la rigidité initale avec la première 
% valeur de rigidité, par défaut: 1
% alpha_p = 0.0; % pourcentage de rigidite tengentielle prp a celle initiale
FY = XE*alpha_e*w_ini^2; % seuil d'elasticite de force 
% comportement pour le modele Mazars
At = alpha_p/alpha_e;
Bt = 0.0;
Ac = alpha_p/alpha_e;
Bc = 0.;
Y0 = XE;
%% la reponse
dt_cal=t_cal(2)-t_cal(1); % pas de temps
u=zeros(1,length(t_cal)); % deplacement
v=zeros(1,length(t_cal)); % vitesse relatif
a=zeros(1,length(t_cal)); % acceleration 
du=zeros(1,length(t_cal)-1); % increment du deplacement
u_p1=zeros(1,length(t_cal)); % deplacement plastique pour oscillateur elastoplastique
%  l'origine 
ETAT=zeros(1,length(t_cal));% Etat de la rigidite du comportement elastique parfaitement plastique
alpha1=zeros(1,length(t_cal)); % alpha de la rigidite du comportement elastique parfaitement plastique
alpha2=zeros(1,length(t_cal)); % alpha de la rigidite du comportement parfaitement endommageant
fint1=zeros(1,length(t_cal)-1);% force elastique parfaitement plastique
fint2=zeros(1,length(t_cal)-1);% force parfaitement endommageant
f_charge=zeros(1,length(t_cal));%  inegalite clausius duhem
Y2 = XE*ones(1,length(t_cal)-1) ; % seuil d'endommagement
D2 = zeros(1,length(t_cal));%  Valeur D indice d'endommagement
%% l'etat initial
u(1,1)= u_0;
v(1,1)= v_0;
alpha2(1,1)= 1;
alpha1(1,1)=1;
ETAT(1,1)= 0;
fint1(1,1) = alpha1(1,1)*w_ini^2*u(1,1);
fint2(1,1) = alpha2(1,1)*w_ini^2*u(1,1);
a(1,1)=-acc_cal(1)-2*xi_ini*w_ini*v(1,1)-ki*fint1(1,1)-(1-ki)*fint2(1,1);
%%
for i=1:(length(t_cal)-1)
    %% schema beta newmark
    [du_predi,u_predi,v_predi,a_predi]=predicteur_betanewmark_20fev(ki,acc_cal(1,i+1),fint1(1,i),fint2(1,i),a(1,i),v(1,i),u(1,i),dt_cal,xi_ini,w_ini,ETAT(1,i),alpha2(1,i));
    du(1,i)= du_predi;
    u(1,i+1)= u_predi;
    v(1,i+1)= v_predi;
    a(1,i+1)= a_predi;
    %% appel de la fonction de correction
    % l'oscillateur elastique parfaitement plastique 
    [fint_co1,alpha_co1,u_p_co1,ETAT_co1,f_charge1,f_charge_co1,signe1] = correcteur_10nov(du(1,i),ETAT(1,i),u_p1(1,i),fint1(1,i),FY,w_ini,alpha_p);
    fint1(1,i+1)=fint_co1;
    alpha1(1,i+1)=alpha_co1;
    u_p1(1,i+1) = u_p_co1;
    ETAT(1,i+1) = ETAT_co1;
    f_charge(1,i+1) = f_charge_co1;
    % l'oscillateur parfaitement endommageant
    [fint_co2,alpha_co2,Y_co2,D_co2]=correcteur_PE_28oct(u(1,i),du(1,i),alpha2(1,i),Y2(1,i),fint2(1,i),w_ini,At,Ac,Bt,Bc,Y0);
    fint2(1,i+1)=fint_co2;
    alpha2(1,i+1)=alpha_co2;
    Y2(1,i+1)=Y_co2;
    D2(1,i+1)=D_co2;
end
% force des 2 ressorts
fint = ki*fint1+(1-ki)*fint2;
alpha = ki*alpha1+(1-ki)*alpha2;
% % figure
% figure,
% subplot(2,2,4),plot(u,fint1)
% legend('elastique parfaitement plastique vs u')
% subplot(2,2,3),plot(u,fint2)
% legend('parfaitement endommageant vs u')
% subplot(2,2,2),plot(u,fint)
% legend('fint vs u')
% subplot(2,2,1),plot(t_cal,alpha)
% legend('ratio de frequence vs time')

