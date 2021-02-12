classdef waveletSave
    %WAVELETSAVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        folder = 'CWT\big memory save\data';
    end
    
    properties
        ID
        Freqs
        MotherWavelet
        Q
    end
    
    methods
        function obj = waveletSave(Freqs, MotherWavelet, Q)
            %WAVELETSAVE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Freqs = Freqs;
            obj.MotherWavelet = MotherWavelet;
            obj.Q = Q;
            
            IDs = waveletSave.getWvltSavesIDs();
            obj.ID = max(IDs) + 1;
            
            obj.createSaveFolder();
        end
        
        function saveCWTfreqK(obj, CWTk, Kfreq)
            CWTsaveFile = [waveletSave.folder, sprintf('\save%d\CWT%d', [obj.ID, Kfreq])];
            save(CWTsaveFile, 'CWTk');
        end
        
        function CWTk = getCWTfreqK(obj, Kfreq)
            CWTsaveFile = [waveletSave.folder, sprintf('\save%d\CWT%d', [obj.ID, Kfreq])];
            load(CWTsaveFile);
        end
    end
    
    methods (Access = private)
        function createSaveFolder(obj)
            saveFolderName = sprintf('save%d', obj.ID);
            if exist([waveletSave.folder, '\', saveFolderName], 'dir')
                error('save file already exists');
            end
            
            mkdir(waveletSave.folder, saveFolderName);
        end
    end
    
    methods (Static, Access = private)
        function IDsList = getWvltSavesIDs()
            foldersList = dir(waveletSave.folder);
            foldersNames = {foldersList.name};
            IDsList = [];
            for k = 1:length(foldersNames)
                folderName = foldersNames{k};
                if length(folderName) >= 4 && strcmp(folderName(1:4), 'save')
                    IDsList(end+1) = str2double(folderName(5:end));
                end
            end
        end
    end
end

