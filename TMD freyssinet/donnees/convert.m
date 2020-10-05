folder = 'TMD freyssinet/donnees/';
filePath = [folder, 'Reponse TMD seul.xlsm'];


%% détection des noms des essais

opts = spreadsheetImportOptions('Sheet', 'Synthesis');
sheets = readtable(filePath,opts);

sheets = table2cell(sheets);
sheets = sheets(2:end);
k = 1;
while k < length(sheets) && ~isempty(sheets{k})
    k = k+1;
end
sheets = sheets(1:k-1);

%% enregistrement de chaque essai

for ks = 1:length(sheets)
    sheetName = sheets{ks};
    
    % lecture des donnes
    sheetTable = xlsread(filePath, sheetName);
    
    ki = 1;
    while ki < size(sheetTable, 1) && isnan(sheetTable(ki, 1))
        ki = ki+1;
    end
    kf = ki;
    while kf < size(sheetTable, 1) && ~isnan(sheetTable(kf, 1))
        kf = kf+1;
    end
    sheetTable = sheetTable(ki:kf-1, 1:end);
    
    t = transpose(sheetTable(1:end, 1));
    X = transpose(sheetTable(1:end, 2:4));
    tableExcel = sheetTable;
    
    essai = sheetName;
    
    % enregistrement
    saveFileName = [num2str(ks), ' - ', sheetName];
    save([folder, saveFileName], 't', 'X', 'essai', 'tableExcel');
end