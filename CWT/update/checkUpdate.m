function upToDate = checkUpdate()
%CHECKUPDATE Summary of this function goes here
%   Detailed explanation goes here

urlProject = 'https://github.com/raphaelcarpine/ondelettes/commits/master/CWT';


%% check every hour
persistent lastTimeCheck; % vide si pas vérifié, sinon [temps_verif, resultat_verif]

if ~isempty(lastTimeCheck) && now <= lastTimeCheck(1) + 1/24 % verification toutes les heures
    upToDate = lastTimeCheck(2);
    return
elseif ~isempty(lastTimeCheck)
    lastTimeCheck = [];
end

%%
try
    % date de derniere MAJ
    pathFolder = mfilename('fullpath');
    pathFolder = fileparts(pathFolder);
    pathCWT = fileparts(pathFolder);
    CWTfolder = dir(pathCWT);
    dateCreation = [];
    for k = 1:length(CWTfolder)
        if strcmp(CWTfolder(k).name, '.')
            dateCreation = CWTfolder(k).datenum;
            break
        end
    end
    if isempty(dateCreation)
        error(' ');
    end
    
    % date de dernier release
    data = webread(urlProject, weboptions('RequestMethod', 'get'));
    dates = strfind(data, 'datetime="');
    if length(dates) < 1
        error(' ');
    end
    dates = dates + length('datetime="');
    dateUpdate = -inf;
    for k = 1:length(dates)
        datesK = data(dates(k):dates(k)+18);
        datesK = strrep(datesK, 'T', ' ');
%         disp(datesK);
        datesK = datenum(datesK);
        if datesK > dateUpdate
            dateUpdate = datesK;
        end
    end
    
    % save
    if dateUpdate > dateCreation
        upToDate = false;
    else
        upToDate = true;
    end
    
    lastTimeCheck = [now, upToDate];
    
catch
    warning('Software update check failed');
    upToDate = true;
end

end

