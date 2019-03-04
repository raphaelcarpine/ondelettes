function[u,v,a]=elastic_27jan(t,ag,w,xi)
gama1=1/2;
beta1=1/4;
t_cal=t;
acc_cal=ag;
dt_cal=t_cal(2)-t_cal(1);
%% caracteristique de l'oscillateur 
w_ini=w; % frequence = 10 Hz
xi_ini=xi;
%% 
u=zeros(length(t_cal),length(w_ini));
v=zeros(length(t_cal),length(w_ini));
a=zeros(length(t_cal),length(w_ini));
% état initial
u(1,1)=0;
v(1,1)=0;
a(1,1)=-acc_cal(1)-2*xi_ini*w_ini*u(1,1)-w_ini^2*u(1,1);
keff=1/beta1/dt_cal^2+2*xi_ini*w_ini*gama1/beta1/dt_cal+w_ini^2;
% intération
for i=1:length(t_cal)-1
    delta=-acc_cal(i+1)+1/beta1/dt_cal^2*u(i,1)+1/beta1/dt_cal*v(i,1)+(1/2/beta1-1)*a(i,1)...
        +2*xi_ini*w_ini*(gama1/beta1/dt_cal*u(i,1)+(gama1/beta1-1)*v(i,1)+(gama1/2/beta1-1)*dt_cal*a(i,1));
    u(i+1,1)=delta/keff;
    a(i+1,1)=1/beta1/dt_cal^2*(u(i+1,1)-u(i,1))-1/beta1/dt_cal*v(i,1)+(1-1/2/beta1)*a(i,1);
    v(i+1,1)=v(i,1)+(1-gama1)*dt_cal*a(i,1)+gama1*dt_cal*a(i+1,1);
end