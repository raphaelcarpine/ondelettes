function dispLinRegTable()
    %% folders
    if exist('D:\simulations elements finis non lin\data', 'dir')
        ridgeFolderPath = 'D:\simulations elements finis non lin\data\data ridges';
    elseif exist('C:\Users\raphael\Documents\resultats simul diff finies', 'dir')
        ridgeFolderPath = 'C:\Users\raphael\Documents\resultats simul diff finies\data ridges';
    else
        error(' ');
    end

    %% param
    fmin = 1;
    fmax = 3;
    QArrVal = [2, 6, 10];
    QArr = cellfun(@num2str, num2cell(QArrVal), 'UniformOutput', false);
    MotherWaveletArr = {'morlet'};
    ridgeContinuityArr = {'none', 'simple', 'reverse', 'double', 'slope3'};
    signalDerivationArr = {'position', 'speed', 'acceleration', 'd acceleration', 'd² acceleration'};
    signalDerivationArrVal = -2:2;
    noiseLevelArrVal = [inf, 20, 15, 10, 5]; % SNR
    noiseLevelArr = cellfun(@num2str, num2cell(noiseLevelArrVal), 'UniformOutput', false);

    Q = QArrVal(1);
    MotherWavelet = MotherWaveletArr{1};
    ridgeContinuity = ridgeContinuityArr{1};
    signalDerivation = signalDerivationArrVal(1);
    noiseLevel = noiseLevelArrVal(1);

    T = getLinRegTable(fmin, fmax, Q, MotherWavelet, ridgeContinuity, signalDerivation, noiseLevel);

    %% figure
    fig = uifigure;
    paramPan = uipanel('Parent', fig, 'Units', 'normalized', 'Position', [0, 0.9 1 0.1]);
    tablePan = uipanel('Parent', fig, 'Units', 'normalized', 'Position', [0, 0, 1, 0.9]);

    tableObj = uitable(tablePan, 'Data', T, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);

    QMenu = uidropdown(paramPan, 'Items', QArr, 'ItemsData', QArrVal, 'ValueChangedFcn', @updateTable,...
        'Position', [5, 5, 100, 22]);
    MotherWaveletMenu = uidropdown(paramPan, 'Items', MotherWaveletArr, 'ValueChangedFcn', @updateTable,...
        'Position', [125, 5, 100, 22]);
    ridgeContinuityMenu = uidropdown(paramPan, 'Items', ridgeContinuityArr, 'ValueChangedFcn', @updateTable,...
        'Position', [245, 5, 100, 22]);
    signalDerivationMenu = uidropdown(paramPan, 'Items', signalDerivationArr, 'ItemsData', signalDerivationArrVal...
        , 'ValueChangedFcn', @updateTable,...
        'Position', [365, 5, 100, 22]);
    noiseLevelMenu = uidropdown(paramPan, 'Items', noiseLevelArr, 'ItemsData', noiseLevelArrVal, 'ValueChangedFcn', @updateTable,...
        'Position', [485, 5, 100, 22]);

    function updateTable(~, ~)
        Q = QMenu.Value;
        MotherWavelet = MotherWaveletMenu.Value;
        ridgeContinuity = ridgeContinuityMenu.Value;
        signalDerivation = signalDerivationMenu.Value;
        noiseLevel = noiseLevelMenu.Value;

        T = getLinRegTable(fmin, fmax, Q, MotherWavelet, ridgeContinuity, signalDerivation, noiseLevel);

        tableObj = uitable(tablePan, 'Data', T, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
%         tableObj = uitable(tablePan, 'Data', T{:,:}, 'ColumnName', T.Properties.VariableNames, ...
%             'RowName', T.Properties.RowNames, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
    end

    %% 

    function T = getLinRegTable(fmin, fmax, Q, MotherWavelet, ridgeContinuity, signalDerivation, noiseLevel)
        signalDerivationStrs = ["_2integration", "_integration", "", "_derivation", "_2derivation"];
        ridgeFolder = sprintf('ridges_fmin%g_fmax%g_Q%g_%s%s_%s', [fmin, fmax, Q,...
            convertCharsToStrings(MotherWavelet), signalDerivationStrs(signalDerivation+3), convertCharsToStrings(ridgeContinuity)]);
        if noiseLevel < inf
            ridgeFolder = [ridgeFolder, '_noise', num2str(noiseLevel)];
        end
        
        try
            load(fullfile(ridgeFolderPath, ridgeFolder, 'linreg.mat'));
            
            % p-val coeff pos
            p_val_diff_6 = 1 - erf(abs(T.coeffPos-T.coeffPos(6))./(sqrt(2)*sqrt(T.coeffPos_err.^2+T.coeffPos_err(6)^2)/1.96));
            T.p_val_diff_6 = p_val_diff_6;
            
            % numeric precision
            for k = 1:size(T, 1)
                T.F0(k) = sprintf("%.3f", T.F0(k));
                T.F0_err(k) = sprintf("%.3f", T.F0_err(k));
                T.coeffAmpl(k) = sprintf("%.2f", T.coeffAmpl(k));
                T.coeffAmpl_err(k) = sprintf("%.2f", T.coeffAmpl_err(k));
                T.coeffPos(k) = sprintf("%.1f", T.coeffPos(k));
                T.coeffPos_err(k) = sprintf("%.1f", T.coeffPos_err(k));
                T.coeffTemp(k) = sprintf("%.4f", T.coeffTemp(k));
                T.coeffTemp_err(k) = sprintf("%.4f", T.coeffTemp_err(k));
            end
        catch
            T = table();
        end
    end

end