clear all

if exist('D:\simulations elements finis non lin\data', 'dir')
    ridgeFolderPath = 'D:\simulations elements finis non lin\data\data ridges';
elseif exist('C:\Users\raphael\Documents\resultats simul diff finies', 'dir')
    ridgeFolderPath = 'C:\Users\raphael\Documents\resultats simul diff finies\data ridges';
else
    error(' ');
end

Ksimul = [1:36, 101, 106];
Qarr = [2 6 10];
fminmaxArr = [1 3];
MotherWaveletArr = {'morlet'};
signalDerivationArr = -2:0;
ridgeContinuityArr = {'none', 'simple', 'reverse', 'double'};
ridgeContinuityArr = {'slope1'};

%%

[initWaitBar, updateWaitBar, closeWaitBar] = getWaitBar(length(Ksimul) * length(Qarr) * size(fminmaxArr, 1) *...
    length(MotherWaveletArr) * length(signalDerivationArr) * length(ridgeContinuityArr),...
    'displayTime', 0, 'windowTitle', 'Computing ridges');
initWaitBar();

% parpool('threads');

for kkf = 1:length(Ksimul)
    kf = Ksimul(kkf);

    % load data file
    try
        filePath = getResultsFile(kf);
        S = load(filePath, 'A', 'L', 'dx', 'T');
        
        % interpolation
        getYtot = @(Y) [zeros(1, size(Y, 2)); Y; zeros(1, size(Y, 2))];
%         S.Y = getYtot(S.Y);
%         S.V = getYtot(S.V);
        S.A = getYtot(S.A);
        pos_capteurs = S.L/2;
%         Ycapt = getYcapt2(S.Y, pos_capteurs, S.dx);
%         Vcapt = getYcapt2(S.V, pos_capteurs, S.dx);
        Acapt = getYcapt2(S.A, pos_capteurs, S.dx);
%         clear Y V A
        
        dataFileExists = true;
    catch
        dataFileExists = false;
    end
    
    for kq = 1:length(Qarr)
        for kfm = 1:size(fminmaxArr, 1)
            for kw = 1:length(MotherWaveletArr)
                for ksd = 1:length(signalDerivationArr)
                    Q = Qarr(kq);
                    fmin = fminmaxArr(kfm, 1);
                    fmax = fminmaxArr(kfm, 2);
                    MotherWavelet = MotherWaveletArr{kw};
                    ct = 0; % edge effect computed later
                    signalDerivation = signalDerivationArr(ksd);
                    
                    % CWT computation
                    freqs = linspace(fmin, fmax, 100);
                    CWT = [];
                    
                    for kc = 1:length(ridgeContinuityArr)
                        ridgeContinuity = ridgeContinuityArr{kc};
                        
                        % check if file doesnt exist
                        if ~dataFileExists
                            updateWaitBar();
                            continue
                        end
                        
                        % file name
                        [~, fileName] = fileparts(filePath);
                        signalDerivationStrs = ["_2integration", "_integration", "", "_derivation", "_2derivation"];
                        ridgeFolder = sprintf('ridges_fmin%g_fmax%g_Q%g_%s%s_%s', [fmin, fmax, Q,...
                            convertCharsToStrings(MotherWavelet), signalDerivationStrs(signalDerivation+3), convertCharsToStrings(ridgeContinuity)]);
                        if isfile(fullfile(ridgeFolderPath, ridgeFolder, [fileName, '.mat']))
                            updateWaitBar();
                            continue
                        end
                        
                        % ridge continuity name correction
                        if length(ridgeContinuity) >= 5 && strcmp(ridgeContinuity(1:5), 'slope')
                            slopeTimeConst = str2double(ridgeContinuity(6:end));
                            ridgeContinuity = 'slope';
                        else
                            slopeTimeConst = nan;
                        end
                        
                        % CWT
                        if isempty(CWT)
                            CWT = WvltComp(S.T, Acapt, freqs, Q, 'MotherWavelet', MotherWavelet,...
                                'DerivationOrder', signalDerivation, 'DisplayWaitBar', false);
                        end
                        ridge = SingleRidgeExtract(S.T, freqs, CWT, MotherWavelet, Q, ct, ridgeContinuity, slopeTimeConst);
                        Fridge = ridge.freq;
                        Aridge = ridge.val;
                        
                        % save
                        flag = mkdir(fullfile(ridgeFolderPath, ridgeFolder));
                        parsave(fullfile(ridgeFolderPath, ridgeFolder, fileName), Fridge, Aridge);
                        disp(fullfile(ridgeFolder, fileName));
                        
                        % waitbar
                        updateWaitBar();
                    end
                end
            end
        end
    end
end
closeWaitBar();
disp('done');

function parsave(fileName, Fridge, Aridge)
    save(fileName, 'Fridge', 'Aridge');
end