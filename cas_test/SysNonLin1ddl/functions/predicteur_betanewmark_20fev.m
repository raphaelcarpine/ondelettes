% le 20 fev 2016
% predicteur pour le comportement mixte
% thermodynamique plus modele mazars



function [du_predi,u_predi,v_predi,a_predi]=predicteur_betanewmark_20fev(ki,acc_n_1,fint_ini1,fint_ini2,a_ini,v_ini,u_ini,dt_cal,xi_ini,w_ini,ETAT_ini,alpha_co)
%% Schema Newmark _HHT alpha
gama1 = 1/2;
beta1 = 1/4;
keff=1/beta1/dt_cal^2+gama1/beta1/dt_cal*2*xi_ini*w_ini+w_ini^2*(ki*(1-ETAT_ini)+(1-ki)*alpha_co);
du_predi=(1/keff)*(-acc_n_1-ki*fint_ini1-(1-ki)*fint_ini2+2*xi_ini*w_ini*((gama1/beta1-1)*v_ini+dt_cal*(gama1/2/beta1-1)*a_ini)+1/beta1/dt_cal*v_ini+(1/2/beta1-1)*a_ini);
u_predi= u_ini + du_predi;
a_predi=du_predi*1/beta1/dt_cal^2-v_ini*1/beta1/dt_cal+a_ini*(1-1/2/beta1);
v_predi=v_ini+a_ini*(1-gama1)*dt_cal+a_predi*gama1*dt_cal;
