function modesMat = groupModShapes(Freqs, Shapes)
%GROUPSHAPESBYMAC Summary of this function goes here
%   Shapes{k_pont}(:, k_freq)
%   Freqs{k_pont}(k_freq)

limitFreq = 0.2; % difference relative, au dessus, modes différents
limitMAC = 0; % en dessous de ce MAC, modes différents
realMAC = true; % prendre le mac des parties réelles

N = length(Freqs);
[M, kMax] = max(cellfun(@length, Freqs));
n = size(Shapes{1}, 1); % nb de ddl

%% repartition fréquences

% répartition initiale
modesMat = nan(M, 1);
freqsMat = nan(M, 1);
for m = 1:length(Freqs{kMax})
    modesMat(m, 1) = m;
    freqsMat(m, 1) = Freqs{kMax}(m);
end

for k0 = 1:N
    % ajout pont k
    modesMat = [modesMat(:, 1:end-1), nan(M, 1), modesMat(:, end)];
    freqsMat = [freqsMat(:, 1:end-1), nan(M, 1), freqsMat(:, end)];
    for m = 1:length(Freqs{k0})
        modesMat(m, k0) = m;
        freqsMat(m, k0) = Freqs{k0}(m);
    end
    
    while true
        % optimisation
        k =  k0;
        for m1 = M-2:-1:1 % déplacement groupes de modes
            score0 = sum(std(freqsMat, 0, 2, 'omitnan') ./ mean(freqsMat, 2, 'omitnan'), 'omitnan');
            
%             k = randi(size(modesMat, 2));
%             m1 = randi(M-1);
            modesMat(m1:end, k) = modesMat([end, m1:end-1], k);
            freqsMat(m1:end, k) = freqsMat([end, m1:end-1], k);
            score = sum(std(freqsMat, 0, 2, 'omitnan') ./ mean(freqsMat, 2, 'omitnan'), 'omitnan');
            if score > score0
                modesMat([end, m1:end-1], k) = modesMat(m1:end, k);
                freqsMat([end, m1:end-1], k) = freqsMat(m1:end, k);
            end
        end
        for a = 1:10000 % inversions modes voisins
            score0 = sum(std(freqsMat, 0, 2, 'omitnan') ./ mean(freqsMat, 2, 'omitnan'), 'omitnan');
            
%             k = randi(size(modesMat, 2));
            m1 = randi(M-1);
            m2 = m1 + 1;
            m1 = randi(M);
            m2 = randi(M);
            modesMat([m1, m2], k) = modesMat([m2, m1], k);
            freqsMat([m1, m2], k) = freqsMat([m2, m1], k);
            score = sum(std(freqsMat, 0, 2, 'omitnan') ./ mean(freqsMat, 2, 'omitnan'), 'omitnan');
            if score > score0
                modesMat([m1, m2], k) = modesMat([m2, m1], k);
                freqsMat([m1, m2], k) = freqsMat([m2, m1], k);
            end
        end
        
        % ajout mode
        % condition fréquences
        ajoutFreq = max(std(freqsMat, 0, 2, 'omitnan') ./ mean(freqsMat, 2, 'omitnan'), [], 'omitnan') > limitFreq;
        % condition deformees
        ajoutMAC = false;
        for m = 1:M
            for k1 = 1:size(modesMat, 2)-1
                if isnan(modesMat(m, k1))
                    continue
                end
                shape1 = Shapes{k1}(:, modesMat(m, k1));
                for k2 = k1+1:size(modesMat, 2)-1
                    if isnan(modesMat(m, k2))
                        continue
                    end
                    shape2 = Shapes{k2}(:, modesMat(m, k2));
                    if realMAC
                        shape1 = real(shape1);
                        shape2 = real(shape2);
                    end
                    MAC = (shape1' * shape2)^2 / ((shape1'*shape1)*(shape2'*shape2));
                    if MAC < limitMAC
                        ajoutMAC = true;
                        break
                    end
                end
            end
        end
        if ~ajoutFreq && ~ajoutMAC
            break
        else
            modesMat = [modesMat; nan(1, size(modesMat, 2))];
            freqsMat = [freqsMat; nan(1, size(freqsMat, 2))];
            M = M + 1;
        end
    end
    
    % tri par freq moyenne croissante
    [~, Isort] = sort(mean(freqsMat, 2, 'omitnan'));
    modesMat = modesMat(Isort, :);
    freqsMat = freqsMat(Isort, :);
end

modesMat = modesMat(:, 1:end-1);


