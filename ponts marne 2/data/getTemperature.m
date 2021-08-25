function [time, temp] = getTemperature(t0)
%GETTEMPERATURE Summary of this function goes here
%   t0: time near sample, datetime nb

%% file read

filePath = 'ponts marne 2\data\temperature.txt';
fid = fopen(filePath);

time = {[]};
temp = {[]};

date0 = '';

tline = fgetl(fid);
while ischar(tline)
    if isempty(tline)
        time{end+1} = [];
        temp{end+1} = [];
    elseif contains(tline, '/')
        date0 = tline;
    else
        C = strsplit(tline, ' ');
        time{end} = [time{end}, datetime([date0, ' ', C{1}], 'InputFormat', 'dd/MM/yyyy HH.mm',...
            'TimeZone', 'Europe/Paris')];
        temp{end} = [temp{end}, str2double(C{2})];
        0;
    end
    
    tline = fgetl(fid);
end
fclose(fid);

kt = 1; % removing empty sets
while kt <= length(time)
    if isempty(time{kt})
        time = time([1:kt-1, kt+1:end]);
        temp = temp([1:kt-1, kt+1:end]);
    else
        kt = kt+1;
    end
end

%% time selection

kt0 = 1; % closest time set
while kt0 <= length(time)-1 && t0 > time{kt0}(end) + (time{kt0+1}(1) - time{kt0}(end))/2
    kt0 = kt0+1;
end

time = time{kt0};
temp = temp{kt0};

if t0 > time(end)
    time = [];
    temp = [];
end


end

