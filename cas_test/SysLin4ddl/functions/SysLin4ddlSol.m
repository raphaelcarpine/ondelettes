function [phi,freq,xi] = SysLin4ddlSol(KSys)

k1 = KSys(1,1);
k2 = KSys(2,1);
k3 = KSys(3,1);
k4 = KSys(4,1);
m1 = KSys(1,2);
m2 = KSys(2,2);
m3 = KSys(3,2);
m4 = KSys(4,2);
c1 = KSys(1,3);
c2 = KSys(2,3);
c3 = KSys(3,3);
c4 = KSys(4,3);
%%
%mi xi'' + xi (ki + ki+1) + xi-1 (-ki-1) + xi+1 (-ki+1) = 0

K =[k1+k2 ,-k2 ,0 ,0;...
    -k2,k2+k3,-k3,0;...
    0,-k3,k3+k4,-k4;...
    0,0,-k4,k4];    

C =[c1+c2 ,-c2 ,0 ,0;...
    -c2,c2+c3,-c3,0;...
    0,-c3,c3+c4,-c4;...
    0,0,-c4,c4];

M = diag([m1,m2,m3,m4]) ;
%%
[phi0,omeg0]=polyeig(K,C,M);
%%

freq=imag(omeg0(1:2:end)/(2*pi));
xi = -real(omeg0(1:2:end)./(freq*2*pi));
phi = (bsxfun(@times,phi0(:,1:2:end),exp(-1i*angle(phi0(1,1:2:end)))));
[~,I] = max(abs(phi),[],1);
phi = bsxfun(@rdivide,phi,[phi(I(1),1),phi(I(2),2),phi(I(3),3),phi(I(4),4)]);