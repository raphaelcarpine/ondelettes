function[fint_co,alpha_co,Y_co,D]=correcteur_PE_28oct(u_before_pre,du_pre,alpha_before_pre,Y_before_pre,fint_before_pre,w_ini,At,Ac,Bt,Bc,Y0)


%%---------------------------------------------%%
% le 28 oct 2015
% corriger selon le modele Mazars
% adapte au comportement PE
%%---------------------------------------------%%
fint_pre = fint_before_pre + du_pre*alpha_before_pre*w_ini^2;
u_pre = u_before_pre+du_pre;
u_equi = sqrt((max(u_pre,0))^2);
r = max(fint_pre,0)/abs(fint_pre);
A = At*(-2*r^2+3*r)+Ac*(2*r^2-3*r+1);
B = r^2*Bt+(1-r^2)*Bc;
% Y = max(Y_before_pre,u_equi); % variable d'endommagement
Y = max(Y_before_pre,abs(u_pre)); % variable d'endommagement
D = 1-(1-A)*Y0/Y-A*exp(-B*(Y-Y0));
fint_co = (1-D)*w_ini^2*u_pre;
alpha_co = (fint_co-fint_before_pre)/du_pre/w_ini^2;
Y_co=Y;




