function  [umax,vmax,amax]=sdfL(td,acc,xi,dtacc,nsteps)

dt=dtacc;
omega=2*pi/td;
alpha=0.25;
delta=0.50;

umax=0.;
vmax=0.;
amax=0.;
uv=zeros(nsteps,1); 
ua=zeros(nsteps,1); 
ut=zeros(nsteps,1); 

% mass=1.;
% Kel=mass*omega*omega;
Kel=1.;
mass=Kel/omega/omega;
c1=2*xi*mass*omega;

kef=Kel+c1*delta/(alpha*dt)+mass/(alpha*dt*dt);

a1=mass/(alpha*dt)+(delta/alpha)*c1;
a2=mass/(alpha*2)+dt*c1*(delta/(2*alpha)-1);

for it=2:nsteps
    
    dacc=acc(it)-acc(it-1);

    Ref=a1*uv(it-1)+a2*ua(it-1)-mass*dacc;

    du=Ref/kef;

    ut(it)=ut(it-1)+du;

    dv=(delta/(dt*alpha))*du-...
        (delta/alpha)*uv(it-1)+dt*ua(it-1)*(1-(delta/(2*alpha)));
    da=...
    (1./(dt*dt*alpha))*du-(1./(alpha*dt))*uv(it-1)-ua(it-1)*(1./(2*alpha));

    uv(it)=uv(it-1)+dv;
    ua(it)=ua(it-1)+da;

    if abs(ut(it)) > umax; umax=abs(ut(it)); end;
    if abs(uv(it)) > vmax; vmax=abs(uv(it)); end;
    if abs(ua(it)+acc(it)) > amax; amax=abs(ua(it)+acc(it)); end;

end;

% hold on
% plot(ua)
% plot(ua+acc,'r-')
