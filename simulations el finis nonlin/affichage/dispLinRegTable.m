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
    QArr = {2, 6, 10};
    MotherWaveletArr = {'morlet'};
    ridgeContinuityArr = {'none', 'simple', 'reverse', 'double'};
    signalDerivationArr = {-2, -1, 0};

    % Q = QArr{1};
    % MotherWavelet = MotherWaveletArr{1};
    % ridgeContinuity = ridgeContinuityArr{1};
    % signalDerivation = signalDerivationArr{1};

    % T = getLinRegTable(fmin, fmax, Q, MotherWavelet, ridgeContinuity, signalDerivation);

    %% figure
    fig = uifigure;
    paramPan = uipanel('Parent', fig, 'Position', [0, 0.9, 1, 0.1]);
    tablePan = uipanel('Parent', fig, 'Position', [0, 0, 1, 0.9]);

    % tableObj = uitable(tablePan, 'Data', T{:,:}, 'ColumnName', T.Properties.VariableNames, ...
    %     'RowName', T.Properties.RowNames, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);

    QMenu = uidropdown(paramPan, 'Items', QArr, 'ValueChangedFcn', @updateTable);
    MotherWaveletMenu = uidropdown(paramPan, 'Items', MotherWaveletArr, 'ValueChangedFcn', @updateTable);
    ridgeContinuityMenu = uidropdown(paramPan, 'Items', ridgeContinuityArr, 'ValueChangedFcn', @updateTable);
    signalDerivationMenu = uidropdown(paramPan, 'Items', signalDerivationArr, 'ValueChangedFcn', @updateTable);

    function updateTable(~, ~)
        Q = QMenu.Value;
        MotherWavelet = MotherWaveletMenu.Value;
        ridgeContinuity = ridgeContinuityMenu.Value;
        signalDerivation = signalDerivationMenu.Value;

        T = getLinRegTable(fmin, fmax, Q, MotherWavelet, ridgeContinuity, signalDerivation);

        tableObj = uitable(tablePan, 'Data', T{:,:}, 'ColumnName', T.Properties.VariableNames, ...
            'RowName', T.Properties.RowNames, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
    end

    %% 

    function T = getLinRegTable(fmin, fmax, Q, MotherWavelet, ridgeContinuity, signalDerivation)
        signalDerivationStrs = ["_2integration", "_integration", "", "_derivation", "_2derivation"];
        ridgeFolder = sprintf('ridges_fmin%g_fmax%g_Q%g_%s%s_%s', [fmin, fmax, Q,...
            convertCharsToStrings(MotherWavelet), signalDerivationStrs(signalDerivation+3), convertCharsToStrings(ridgeContinuity)]);
        
        load(fullfile(ridgeFolderPath, ridgeFolder, 'linreg.mat'));
    end

end