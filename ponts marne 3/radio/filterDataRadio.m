folder = 'C:\Users\carpine\Documents\projets\ponts marne\reprise operations 2022\donnees\radio';
load(fullfile(folder, 'esbly0107_aprem.mat'));

%% reechantillonage

% 128.0000
% 200.0058
% 
X = resample(X.', 128, 100).';
T = (0:size(X, 2)-1)/128;

%% passages devant le laser

% X(:, T>5490 & T<5493) = nan;
% X(:, T>1776 & T<1785) = nan;
% 
% % decoupage
% k1i = 1;
% k1f = 1;
% while ~isnan(X(1, k1f+1))
%     k1f = k1f + 1;
% end
% 
% k2i = k1f+1;
% while isnan(X(1, k2i))
%     k2i = k2i + 1;
% end
% k2f = k2i;
% while ~isnan(X(1, k2f+1))
%     k2f = k2f + 1;
% end
% 
% k3i = k2f+1;
% while isnan(X(1, k3i))
%     k3i = k3i + 1;
% end
% k3f = length(T);

%% filtrage

Tmean = 30;
for kl = 1:size(X, 1)
%     X(kl, k1i:k1f) = X(kl, k1i:k1f) - getSmoothSignal(T(k1i:k1f), X(kl, k1i:k1f), 'gaussian', Tmean);
%     X(kl, k2i:k2f) = X(kl, k2i:k2f) - getSmoothSignal(T(k2i:k2f), X(kl, k2i:k2f), 'gaussian', Tmean);
%     X(kl, k3i:k3f) = X(kl, k3i:k3f) - getSmoothSignal(T(k3i:k3f), X(kl, k3i:k3f), 'gaussian', Tmean);
    
    X(kl, :) = X(kl, :) - getSmoothSignal(T, X(kl, :), 'gaussian', Tmean);
end


%% enregistrement

Tradio = T;
Xradio = X;

save(fullfile(folder, 'esbly0107_aprem_resampled_filtered.mat'), 'Tradio', 'Xradio', 'chDist');









