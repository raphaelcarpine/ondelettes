ridgeFolderPath = 'C:\Users\carpine\Documents\projets\simulations elements finis non lin\data\data ridges';

Nsimul = 26;
Qarr = [2 6 10];
fminmaxArr = [1 3];
MotherWaveletArr = {'morlet'};
ridgeContinuityArr = [0 1];
signalDerivation = -2;

%%

[initWaitBar, updateWaitBar, closeWaitBar] = getWaitBar(Nsimul * length(Qarr) * size(fminmaxArr, 1) * length(MotherWaveletArr),...
    'displayTime', 0, 'windowTitle', 'Computing ridges');
initWaitBar();

for kf = 101%1:Nsimul
    % load data file
    try
        filePath = getResultsFile(kf);
        load(filePath);
        
        % interpolation
        getYtot = @(Y) [zeros(1, size(Y, 2)); Y; zeros(1, size(Y, 2))];
        Y = getYtot(Y);
        V = getYtot(V);
        A = getYtot(A);
        pos_capteurs = L/2;
        Ycapt = getYcapt2(Y, pos_capteurs, dx);
        Vcapt = getYcapt2(V, pos_capteurs, dx);
        Acapt = getYcapt2(A, pos_capteurs, dx);
        clear Y V A
        
        dataFileExists = true;
    catch
        dataFileExists = false;
    end
    
    for kq = 1:length(Qarr)
        for kfm = 1:size(fminmaxArr, 1)
            for kw = 1:length(MotherWaveletArr)
                Q = Qarr(kq);
                fmin = fminmaxArr(kfm, 1);
                fmax = fminmaxArr(kfm, 2);
                MotherWavelet = MotherWaveletArr{kw};
                
                % check if file already exists
                if ~dataFileExists
                    updateWaitBar();
                    continue
                end
                [~, fileName] = fileparts(filePath);
                signalDerivationName = ["_2integration", "_integration", "", "_derivation", "_2derivation"];
                ridgeFolder = sprintf('ridges_fmin%g_fmax%g_Q%g_%s%s', [fmin, fmax, Q,...
                    convertCharsToStrings(MotherWavelet), signalDerivationName(signalDerivation+3)]);
                if isfile(fullfile(ridgeFolderPath, ridgeFolder, [fileName, '.mat']))
                    updateWaitBar();
                    continue
                end
                
                % CWT computation
                freqs = linspace(fmin, fmax, 100);
                CWT = WvltComp(T, Acapt, freqs, Q, 'MotherWavelet', MotherWavelet,...
                    'DerivationOrder', signalDerivation, 'DisplayWaitBar', false);
                freqs = [nan, freqs, nan];
                CWT = [zeros(1, size(CWT, 2)); CWT; zeros(1, size(CWT, 2))];
                
                for kc = 1:length(ridgeContinuityArr)
                    ridgeContinuity = ridgeContinuityArr(kc);
                    
                    % ridge computation
                    if ridgeContinuity
                        Fridge = nan(1, size(CWT, 2));
                        [~, Fridge(1)] = max(abs(CWT(:, 1)));
                        localMax = abs(CWT(1:end-1, :)) < abs(CWT(2:end, :));
                        localMax = localMax(1:end-1, :) & ~localMax(2:end, :);
                        for kt = 2:size(CWT, 2)
                            localMaxFreq = find(localMax(:, kt)) + 1;
                            [~, closestLocalMax] = min(abs(localMaxFreq - Fridge(kt-1)));
                            Fridge(kt) = localMaxFreq(closestLocalMax);
                        end
                    else
                        [~, Fridge] = max(abs(CWT), [], 1);
                    end
                    [Fridge, Aridge] = localMax3Points(freqs([Fridge-1; Fridge; Fridge+1]),...
                        CWT([Fridge-1; Fridge; Fridge+1] + [1;1;1] * (0:size(CWT, 2)-1)*size(CWT, 1)));
                    
                    % save
                    if ridgeContinuity
                        ridgeFolder2 = [ridgeFolder, '_continuous'];
                    else
                        ridgeFolder2 = ridgeFolder;
                    end
                    if ~exist(fullfile(ridgeFolderPath, ridgeFolder2), 'dir')
                        mkdir(fullfile(ridgeFolderPath, ridgeFolder2));
                    end
                    save(fullfile(ridgeFolderPath, ridgeFolder2, fileName), 'Fridge', 'Aridge');
                end
                
                updateWaitBar();
            end
        end
    end
end
closeWaitBar();