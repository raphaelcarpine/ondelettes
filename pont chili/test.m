clear all

%% data
Fs = 100;
DOFs = 1:9;
load('pont chili/data/Event_2018-05-02_13-10-35_969');
% load('pont chili/data/Periodic_2018-04-26_11-59-54_969');
% load('pont chili/data/Periodic_2018-05-10_17-59-55_969');
removeMean = true;
autocorrelation = true;

dt = 1/Fs;
t = dt * (0:size(X, 2)-1);

X = X(DOFs, :);
X = X*9.81;


%%

% if removeMean
%     for dof = 1:size(X, 1)
%         X(dof, :) = X(dof, :) - mean(X(dof, :));
%     end
% end

% if autocorrelation
%     Rx = nan(size(X));
%     for dof = 1:size(X, 1)
%         RxDof = xcorr(X(dof, :), 'biased') / var(X(dof, :));
%         Rx(dof, :) = RxDof(ceil(length(RxDof)/2):end);
%     end
% end

%% plots

fig = figure;
ax = axes(fig);
plts = plot(t, X);
xlabel(ax, 'Time [s]');
ylabel(ax, 'Acceleration [m/s²]');

% if autocorrelation
%     fig = figure;
%     ax = axes(fig);
%     pltsAutocorr = plot(t, Rx);
%     xlabel(ax, 'Time [s]');
%     ylabel(ax, 'Ra');
% end


%% wavelet

Q = 10;
NbMaxRidges = 1;
NbMaxParallelRidges = 1;
fmin = 1;
fmax = 10;

WaveletMenu('WaveletPlot', plts, 'Q', Q, 'fmin', fmin, 'fmax', fmax, 'AutocorrelationMode', autocorrelation,...
    'MaxRidges', NbMaxRidges, 'MaxParallelRidges', NbMaxParallelRidges, 'RemoveMean', removeMean);


% if autocorrelation
%     WaveletMenu('WaveletPlot', pltsAutocorr, 'Q', Q, 'fmin', fmin, 'fmax', fmax, 'MultiSignalMode', true,...
%     'MaxRidges', NbMaxRidges, 'MaxParallelRidges', NbMaxParallelRidges);
% end

















% 
% 
% %%
% 
% %(K-w2*M)*phi = 0 => (M\K-w2*Id)*phi = 0 => M\K*phi = w2*phi
% 
% mu = 1;
% deltaX = 1;
% N = 10;
% masse_essieu = 10;
% 
% 
% M = eye(N)*mu*deltaX;
% 
% 
% K = 6*eye(N);
% for index = 1:N-1
%     K(index, index+1) = -4;
%     K(index+1 , index) = -4;
% end
% for index = 1:N-2
%     K(index, index+2) = 1;
%     K(index+2 , index) = 1;
% end
% 
% 
% 
% xessieu = 5;
% DDLessieu = round(xessieu/deltaX);
% Mtilde = M;
% Mtilde(DDLessieu, DDLessieu) = Mtilde(DDLessieu, DDLessieu) + masse_essieu;
% 
% 
% 
% [Vtilde, Dtilde] = eig(Mtilde\K);
% 
% for mode = 1:N
%     deformee = Vtilde(:, mode);
%     m_modal = traspose(deformee)*Mtilde*deformee;
%     k_modal = traspose(deformee)*K*deformee;
%     c_modal = traspose(deformee)*C*deformee;
%     
%     
% 











































