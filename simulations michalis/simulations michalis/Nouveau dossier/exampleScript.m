clear all; clc; close all;


% READ GROUND MOTION RECORDS

recsPath = fullfile(fileparts(matlab.desktop.editor.getActiveFilename), 'recs');

% [acc,vel]=ReadRecord('184057.AT2',recsPath,'NGA');
% plot((1:acc.nsteps)*acc.dtacc,acc.rec)
% 
% 
% [acc,vel]=ReadRecord('PL_1_SN.acc',recsPath,'NGA2');
% plot((1:acc.nsteps)*acc.dtacc,acc.rec)
% 
% 
% [acc,vel]=ReadRecord('NGA_no_28_C12050.AT2',recsPath,'NGA');
% plot((1:acc.nsteps)*acc.dtacc,acc.rec)
% 
% 
% [acc,vel]=ReadRecord('AssisiStallone_X.txt',recsPath,'plain3');
% plot((1:acc.nsteps)*acc.dtacc,acc.rec)

[acc,vel]=ReadRecord('sv23.th',recsPath,'Silva');
plot((1:acc.nsteps)*acc.dtacc,acc.rec)


% GET SPECTRAL VALUES FOR A SINGLE PERIOD VALUE (T=1sec)
[umax,vmax,amax]=sdfL(1,acc.rec,0.05,acc.dtacc,acc.nsteps)


% RESPONSE SPECTRUM (damping = 5%)
[td,amax,vmax,umax] = spectrum(acc.rec,acc.dtacc,0.05);
plot(td,amax)


% RUN A NONLINEAR SDOF
Tn = 1; %sec
xi = 0.05;
dt = acc.dtacc;
hardening = 0.01; % post yield over elastic modulus
fy = 0.01;
nsteps = acc.nsteps;
niter  = 10;
sfactor = 1; % record scaling factor
iplot = 1; % draw hysteretic plot

[umax,vmax,amaxRel,amaxTot,fsmax,u,fs]=...
    sdfNL(Tn,xi,dt,hardening,fy,acc.rec,nsteps,niter,sfactor,iplot);