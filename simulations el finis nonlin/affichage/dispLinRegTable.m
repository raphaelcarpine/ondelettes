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
    ridgeContinuityArr = {'none', 'simple', 'reverse', 'double'};
    signalDerivationArr = {'position', 'speed', 'acceleration'};
    signalDerivationArrVal = [-2, -1, 0];

    Q = QArrVal(1);
    MotherWavelet = MotherWaveletArr{1};
    ridgeContinuity = ridgeContinuityArr{1};
    signalDerivation = signalDerivationArrVal(1);

    T = getLinRegTable(fmin, fmax, Q, MotherWavelet, ridgeContinuity, signalDerivation);

    %% figure
    fig = uifigure;
    paramPan = uipanel('Parent', fig, 'Units', 'normalized', 'Position', [0, 0.9 1 0.1]);
    tablePan = uipanel('Parent', fig, 'Units', 'normalized', 'Position', [0, 0, 1, 0.9]);

    tableObj = uitable(tablePan, 'Data', T{:,:}, 'ColumnName', T.Properties.VariableNames, ...
        'RowName', T.Properties.RowNames, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);

    QMenu = uidropdown(paramPan, 'Items', QArr, 'ItemsData', QArrVal, 'ValueChangedFcn', @updateTable,...
        'Position', [5, 5, 100, 22]);
    MotherWaveletMenu = uidropdown(paramPan, 'Items', MotherWaveletArr, 'ValueChangedFcn', @updateTable,...
        'Position', [125, 5, 100, 22]);
    ridgeContinuityMenu = uidropdown(paramPan, 'Items', ridgeContinuityArr, 'ValueChangedFcn', @updateTable,...
        'Position', [245, 5, 100, 22]);
    signalDerivationMenu = uidropdown(paramPan, 'Items', signalDerivationArr, 'ItemsData', signalDerivationArrVal...
        , 'ValueChangedFcn', @updateTable,...
        'Position', [365, 5, 100, 22]);

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