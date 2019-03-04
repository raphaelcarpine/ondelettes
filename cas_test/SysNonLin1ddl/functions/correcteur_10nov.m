function[fint_co,alpha_co,u_p_co,ETAT_co,f_charge,f_charge_co,signe]=correcteur_10nov(du_pre,ETAT_before_pre,u_p_bef_pre,fint_before_pre,FY,w_ini,alpha_p)
% 29 oct 2015
% elasto plastique avec ecrouissage
% thermodynamique
% corriger comme la correction faite pour le programme utilisant le seuil
% evolue
% le resultat pour ce programme est satisfait
dup_pre = du_pre*ETAT_before_pre;
up_pre = u_p_bef_pre+dup_pre;
fint_pre = fint_before_pre+(du_pre-dup_pre)*w_ini^2;
f_charge = abs(fint_pre-(alpha_p/(1-alpha_p))*w_ini^2*up_pre)-FY;
if ETAT_before_pre == 0 % elastique
    if f_charge < 0 % continuer en elastique
        dup_co = dup_pre;
        fint_co = fint_pre;
        u_p_co = u_p_bef_pre+dup_co;
        f_charge_co = abs(fint_co-(alpha_p/(1-alpha_p))*w_ini^2*u_p_co)-FY;
        if f_charge_co <0
            ETAT_co = 0;
            alpha_co = 1;
        else
            ETAT_co = 1-alpha_p;
            alpha_co = alpha_p;
        end
        
    else % en train d'aller a la plastique
        if sign(fint_pre-(alpha_p/(1-alpha_p))*w_ini^2*up_pre) >= 0 
            F_limite = FY +(alpha_p/(1-alpha_p))*w_ini^2*u_p_bef_pre;
        else
            F_limite = -FY +(alpha_p/(1-alpha_p))*w_ini^2*u_p_bef_pre;
        end
        delta_ue = ((F_limite-fint_before_pre)*du_pre)/((du_pre-dup_pre)*w_ini^2);
        delta_up = du_pre-delta_ue;
        dup_co= (1-alpha_p)*delta_up;
        u_p_co = u_p_bef_pre+dup_co;
        if sign(fint_pre-(alpha_p/(1-alpha_p))*w_ini^2*up_pre) >= 0 
            fint_co = FY +(alpha_p/(1-alpha_p))*w_ini^2*u_p_co;
        else
            fint_co = -FY +(alpha_p/(1-alpha_p))*w_ini^2*u_p_co;
        end
        ETAT_co = 1- alpha_p;
        alpha_co = alpha_p;
        f_charge_co = abs(fint_co-(alpha_p/(1-alpha_p))*w_ini^2*u_p_co)-FY;
        
    end
else % en plastique
    if sign(fint_pre-(alpha_p/(1-alpha_p))*w_ini^2*up_pre)*dup_pre >=0 % continuer a plastifier
        dup_co = (1-alpha_p)*du_pre;
        u_p_co = u_p_bef_pre+dup_co;
        if sign(fint_pre-(alpha_p/(1-alpha_p))*w_ini^2*up_pre) >= 0 
            fint_co = FY +(alpha_p/(1-alpha_p))*w_ini^2*u_p_co;
        else
            fint_co = -FY +(alpha_p/(1-alpha_p))*w_ini^2*u_p_co;
        end
        ETAT_co = 1- alpha_p;
        alpha_co = alpha_p;
        f_charge_co = abs(fint_co-(alpha_p/(1-alpha_p))*w_ini^2*u_p_co)-FY;
   else % decharge depuis la plasticite
        dup_co = 0 ; %0
        fint_co = fint_before_pre+ dup_pre*w_ini^2;
        u_p_co = u_p_bef_pre+dup_co;
        f_charge_co = abs(fint_co-(alpha_p/(1-alpha_p))*w_ini^2*u_p_co)-FY;
        if f_charge_co <0
            ETAT_co = 0;
            alpha_co = 1;
        else 
            ETAT_co = 1-alpha_p;
            alpha_co =alpha_p;
        end
   end
   
    
end
f_charge_co = abs(fint_co-(alpha_p/(1-alpha_p))*w_ini^2*u_p_co)-FY;
signe=sign(fint_pre-(alpha_p/(1-alpha_p))*w_ini^2*up_pre)*dup_pre;