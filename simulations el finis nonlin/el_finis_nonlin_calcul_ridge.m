ridgeFolderPath = 'C:\Users\carpine\Documents\projets\simulations elements finis non lin\data\data ridges';

Nsimul = 26;
Qarr = [2 6 10];
fminmaxArr = [1 3];
MotherWaveletArr = {'morlet'};
ridgeContinuityArr = {'none', 'simple'};
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
                ct = 0; % edge effect computed later
                
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
                
                for kc = 1:length(ridgeContinuityArr)
                    ridgeContinuity = ridgeContinuityArr{kc};
                    
                    if length(ridgeContinuity) >= 5 && strcmp(ridgeContinuity(1:5), 'slope')
                        slopeTimeConst = str2double(ridgeContinuity(6:end));
                        ridgeContinuity = 'slope';
                    else
                        slopeTimeConst = nan;
                    end
                    
                    ridge = SingleRidgeExtract(T, freqs, CWT, MotherWavelet, Q, ct, ridgeContinuity, slopeTimeConst);
                    Fridge = ridge.freq;
                    Aridge = ridge.val;
                    
                    % save
                    ridgeFolder2 = [ridgeFolder, '_', ridgeContinuity];
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