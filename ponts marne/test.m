clear all

%% data
bridge = 'trilbardou';
Fs = 50;
load('ponts marne/data/trilbardou1506/ambiant4_10min50Hz26C');
meanToZero = true;
correctClipping = true;

dt = 1/Fs;
t = dt * (0:size(X, 2)-1);


% X = X([2], :);

%% sensors bugs

boundsX = [-20, 0];

if correctClipping
    clippingChanels = [];
    for dof = 1:size(X, 1)
        clippingInds = X(dof, :) < boundsX(1) | X(dof, :) > boundsX(2);
        if any(clippingInds)
            clippingChanels = [clippingChanels, dof];
            X(dof, clippingInds) = nan;
        end
    end
    
    for dof = clippingChanels
        for kt = 2:size(X, 2)
            if isnan(X(dof, kt))
                X(dof, kt) = X(dof, kt-1);
            end
        end
        for kt = size(X, 2)-1:-1:1
            if isnan(X(dof, kt))
                X(dof, kt) = X(dof, kt+1);
            end
        end
    end
    
    if ~isempty(clippingChanels)
        warning(['clipping in channels ', num2str(clippingChanels)]);
    end
end

%%

if meanToZero
    for dof = 1:size(X, 1)
        X(dof, :) = X(dof, :) - mean(X(dof, :));
    end
end


%% plots

fig = figure;
ax = axes(fig);
plts = plot(t, X);
xlabel(ax, 't');
ylabel(ax, 'a');


%% shapes plot

dimensionsShapes

%% wavelet

Q = 10;
NbMaxRidges = 1;
NbMaxParallelRidges = 1;
fmin = 0.5;
fmax = 5;

WaveletMenu('WaveletPlot', plts, 'Q', Q, 'fmin', fmin, 'fmax', fmax, 'AutocorrelationMode', true, 'MaxLagCorr', 100,...
    'MaxRidges', NbMaxRidges, 'MaxParallelRidges', NbMaxParallelRidges, 'RealShapePlot', shapePlotBridge);

















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











































