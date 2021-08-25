%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RUN SDOF OSCILLATOR
%
% created by Michalis Fragiadakis, Dec 2013
% mfrag@mail.ntua.gr
%
% please report any bugs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
includeFile

% RUN A NONLINEAR SDOF
Tn = 1; %sec
xi = 0.05;
dt = acc.dtacc;
hardening = 0.01; % post yield over elastic modulus

nsteps = acc.nsteps;
niter  = 10;
sfactor = 1; % record scaling factor
iplot = 0; % draw hysteretic plot

% CHANGE THIS TO MOVE FROM ELASTIC TO ELASTOPLASTIC
% fy = 0.01;
fy = 0.03;
% fy = 100;

[umax,vmax,amaxRel,amaxTot,fsmax,u,fs,a]=...
    sdfNL(Tn,xi,dt,hardening,fy,acc.rec,nsteps,niter,sfactor,iplot);

% plot hysteretisis
figure()
hold on; grid on; box on;
plot(u,fs,'k.-')

xlabel('displacement','FontSize',18)
ylabel('force','FontSize',18)

t = (1:(nsteps-1))*dt;


% plot acceleration over time
figure()
hold on; grid on; box on;
plt = plot(t,a,'k-');

xlabel('time (sec)','FontSize',18)
ylabel('acceleration','FontSize',18)

% plot force over time
figure()
hold on; grid on; box on;
plot(t,fs,'k-');
xlabel('time (sec)','FontSize',18)
ylabel('force','FontSize',18)

% plot exciattion
figure;
hold on; grid on; box on;
plt2 = plot(dt*(0:length(acc.rec)-1), acc.rec, 'k-');
xlabel('time (sec)','FontSize',18)
ylabel('input','FontSize',18)

% %%
% % Scalogram of acceleration response signal
% Fs = linspace(0.1,5);
% [WvltOut, ctZeroPadding] = WvltComp(t,a, Fs, 5);
% 
% figure()
% [X,Y] = meshgrid(Fs,t);
% Z = abs(WvltOut);
% %contour(X,Y,Z')
% h = pcolor(Y,X,Z');
% set(h, 'EdgeColor', 'none')

% %%
% % Scalogram of displacement response signal
% Fs = linspace(0.1,5);
% [WvltOut, ctZeroPadding] = WvltComp(t,u, Fs, 5);
% 
% figure()
% [X,Y] = meshgrid(Fs,t);
% Z = abs(WvltOut);
% %contour(X,Y,Z')
% h = pcolor(Y,X,Z')
% set(h,'EdgeColor','none')

%% wavelet menu

fmin = 0.1;
fmax = 10;
Q = 5;
WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q);
WaveletMenu('WaveletPlot', plt2, 'fmin', fmin, 'fmax', fmax, 'Q', Q);



