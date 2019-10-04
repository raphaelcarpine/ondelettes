function [t, X] = getData(P, transient)
%GETDATA function to recover data from the data folder
%   P : load index. 0, 6 or 7.
%   transient : transient index. 0 for whole signal, 1, 2, 3... for
%   transients
%   t : time, 1*n matrix
%   X : signal data, 9*n matrix (nine sensors)

P = max(P-4, 1); %0->1 ; 6->2 ; 7->3
transient = transient + 1;

%%
if transient == 1
    load('mur silvia/data/allData/mData.mat');
else
    load('mur silvia/data/allTransients/mData.mat');
end

%%
varNamesPrefix = cell(3, 4); %before the sensor index
varNamesSuffix = cell(3, 4); %after the sensor index

varNamesPrefix{1,1} = 'ch';
varNamesSuffix{1,1} = '_total';

varNamesPrefix{1,2} = 'transient_1_ch';
varNamesSuffix{1,2} = '';

varNamesPrefix{1,3} = 'transient_2_ch';
varNamesSuffix{1,3} = '';

varNamesPrefix{1,4} = 'transient_3_ch';
varNamesSuffix{1,4} = '';

varNamesPrefix{2,1} = 'P6_ch';
varNamesSuffix{2,1} = '_total';

varNamesPrefix{2,2} = 'transient_1_P6_ch';
varNamesSuffix{2,2} = '';

varNamesPrefix{2,3} = 'transient_2_P6_ch';
varNamesSuffix{2,3} = '';

varNamesPrefix{3,1} = 'P7_ch';
varNamesSuffix{3,1} = '_total';

varNamesPrefix{3,2} = 'transient_1_P7_ch';
varNamesSuffix{3,2} = '';

varNamesPrefix{3,3} = 'transient_2_P7_ch';
varNamesSuffix{3,3} = '';

varNamesPrefix{3,4} = 'transient_3_P7_ch';
varNamesSuffix{3,4} = '';

%%
X = cell(1, 9);
for indS = 1:9
    eval(['X{indS} = ', varNamesPrefix{P, transient}, num2str(indS), varNamesSuffix{P, transient}, ';']);
end

%%
X = transpose([X{:}]);
t = deltat * (0:(size(X, 2)-1));



end

