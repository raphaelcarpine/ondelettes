function [td,amax,vmax,umax] = spectrum(aa,dt,damp,td)
% [td,amax,vmax,umax] = spectrum(aa,dt,damp,td)
% Spectrum - calculate record response spectrum
%IN:
% aa    = acceleration
% dt    = time stem
% damp  = damping, xi=0.05 (not 5)
% td (optional) = period values, given as a vector
%
% OUT
% umax, vmax, amax = displacement, velocity, acceleration

if nargin < 4
    td = (1:500)*0.01;
end;


npts=length(aa);

%   Detailed explanation goes here

for j=1:length(td)
    [umax(j),vmax(j),amax(j)]=sdfL(td(j),aa,damp,dt,npts); % new signal
end;