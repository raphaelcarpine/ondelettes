% function de l'intégration par beta newmark
function[v,u]=beta_newmark(t,a,v_ini,u_ini)
v = zeros(length(a),1);
u = zeros(length(a),1);
v(1) = v_ini;
u(1) = u_ini;
beta = 1/4;
gama = 1/2;
dt = t(2) - t(1);
for i=1:length(t)-1
    v(i+1) = v(i) + (1-gama)*dt*a(i)+gama*dt*a(i+1);
    u(i+1) = u(i) + dt*v(i) + (1-2*beta)/2*dt^2*a(i)+beta*dt^2*a(i+1);
end
