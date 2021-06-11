%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BILINEAR HYSTERETIC LAW
%
% created by Michalis Fragiadakis, Dec 2013
% mfrag@mail.ntua.gr
%
% please report any bugs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigs,Es0]=ForceDelta(Fsy,Eo,b,epssP,sigsP,depss)
%
% Eo: young modulus
% b:        strain hardening ratio Et/E
% epssP:    previous strain
% sigsP:    previous stress
% depss:    difference in strain

%..	calculate other fixed material properties
Esh  = b*Eo;
epsy = Fsy/Eo;

%..	calculate current strain
epss   = epssP + depss;

%..	bilinear model
c1   = Esh*epss;
c2   = (Fsy)*(1.-b);
c3   = (Fsy)*(1.-b);
c    = sigsP + Eo*depss;
sigs = max((c1-c2),min((c1+c3),c));
Es0=Esh;
if (abs(sigs-c) < 1.d-10) Es0=Eo; end;

epssP=epss;
sigsP=sigs;