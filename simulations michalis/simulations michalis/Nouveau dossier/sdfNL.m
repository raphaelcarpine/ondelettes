function  [umax,vmax,amaxRel,amaxTot,fsmax,u,fs,a]=...
    sdfNL(Tn,xi,dt,hardening,fy,ff,nsteps,niter,sfactor,iplot)

if nargin<10
    iplot=0;
end
currentdir=cd;

% TEST LINES
% clear
% dt=0.01;
% fy=7.5;
% Kel=10;
% mass=0.2533;
% damp=0.1592;
% totalTime=10;
% ff=forceSine(totalTime,dt,10,0.6); % forceSine(tt,dt,Ap,td2)
% nsteps=totalTime/dt;

% DO SOME PRELIMINARY CALCULATIONS
omega=2*pi/Tn;
mass=1.;
damp=2*xi*mass*omega;
Kel=omega*omega/mass;
ff=sfactor*ff;

alpha=4*mass/dt+2*damp;
beta=2*mass;

k0=Kel; k(1)=Kel;
p(1)=0; v(1)=0; u(1)=0;
a(1)=(p(1)-damp*v(1)-k(1)*u(1))/mass;
dphatP=0;
fs(1)=0;
alpha0=0;
k(2)=k(1);
amax=0;


for t=2:nsteps-1
    
    %     waitbar(t/nsteps)
    dp=ff(t)-ff(t-1);
    dphat=dp+alpha*v(t-1)+beta*a(t-1);
    
    if dphat*dphatP < 0
        if k(t)==hardening*k(t)
            k(t)=k0;
        elseif k(t)==k0
            k(t)=hardening*k(t);
        end
    end
    
    keff=k(t)+2*damp/dt+4*mass/dt/dt;
    
    uTemp=0;
    fsITp=fs(t-1);
    uP=u(t-1);
    uITtot=0;
    
    for iter=1:niter
        
        du=dphat/keff;
        uITtot=uITtot+du; % iterative displacement
        uTemp=uP+uITtot;  % total displacement
        
        [fsIT,k(t+1)]=ForceDelta(fy,k0,hardening,uTemp-du,fsITp,du);
        
        df=fsIT-fsITp+alpha*du/dt;
        dr=dphat-df;
        
        err=dr;
        if err < 0.000001
            break;
        end;
        
        dphat=dr;
        fsITp=fsIT;
    end;
    
    fs(t)=fsIT;
    u(t)=uTemp; %u(t)=u(t-1)+du
    dv=2*(u(t)-u(t-1))/dt-2*v(t-1);
    v(t)=v(t-1)+dv;
    a(t)=(1/mass)*(ff(t)-damp*v(t)-fs(t));
    dphatP=dphat;
    %     if abs(a(t)-ff(t)) > a(t); a(t)=abs(a(t)-ff(t)); end;
end

umax=max(abs(u));
vmax=max(abs(v));
amaxRel=max(abs(a));  % max relative acceleration
fsmax=max(abs(fs));

[i1,i2]=size(ff);
if i1>i2
    ff=ff';
end
ll=min(length(a),length(ff));
amaxTot=max(abs(a(1:ll)-ff(1:ll)));  % max relative acceleration

% save ('plotNL1','u','fs')

if iplot==1
    clf; hold on;
    grid on; box on;
    plot(u,fs,'k--.');
    % comet(u,fs);
end

cd(currentdir);