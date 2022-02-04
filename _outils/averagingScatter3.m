function [X, Y, stdY, K, Xlims] = averagingScatter3(x, y, dx, xIncremenScale, t, DeltaT)
%AVERAGINGSCATTER Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    t = nan;
end

if ~iscolumn(x) && ~isrow(x)
    error('x matrix');
elseif ~iscolumn(y) && ~isrow(y)
    error('y matrix');
elseif length(x) ~= length(y)
    error('length(x) ~= length(y)');
end

%% nan removal

indexes_not_nan = ~isnan(x) & ~isnan(y);
x = x(indexes_not_nan);
y = y(indexes_not_nan);
t = t(indexes_not_nan);


%% time separation

if ~isnan(t)
    t = t - t(1);
    Kt = 0;
    Tindex = nan(size(x));
    for kt = 1:length(t)
        if t(kt) >= length(Kt)*DeltaT
            Kt(end+1) = 0;
        end
        Kt(end) = Kt(end) + 1;
        Tindex(kt) = length(Kt);
    end
%     Nt = floor((t(end)-t(1))/DeltaT); % nb of time intervals
%     
%     kt = floor(length(t)/Nt);
%     rt = mod(length(t), kt);
%     randomKt = 1:Nt;
%     for i = 1:(Nt-rt)
%         randomKt(randi(length(randomKt))) = [];
%     end
%     Kt = kt * ones(1, Nt); % nb of points in time interval
%     Kt(randomKt) = kt+1;
%     
%     Tindex = nan(size(x)); % time interval index
%     for it = 1:length(Kt)
%         Tindex(sum(Kt(1:it-1))+1:sum(Kt(1:it))) = it;
%     end
else
    Kt = ones(size(x));
    Tindex = 1:length(x);
end


%% x separation

switch xIncremenScale
    case 'lin'
    case 'log'
        x = log(x);
    otherwise
        error('');
end

n = floor(x/dx);
ni = min(n);
nf = max(n);
N = nf - ni + 1;

X = ((ni:nf)+1/2) * dx;
Xlims = (ni:nf+1) * dx;
Y0 = cell(N, length(Kt));

for kx = 1:length(x)
    n = floor(x(kx)/dx) - ni + 1;
    Y0{n, Tindex(kx)}(end+1) = y(kx);
end
Y00 = cell(1, N);
for n = 1:N
    for kt = 1:length(Kt)
        Y00{n}(end+1) = mean(Y0{n, kt});
    end
    Y00{n} = Y00{n}(~isnan(Y00{n}));
end

Y = nan(1, N);
K = nan(1, N);
stdY = nan(1, N);

for n = 1:N
    K(n) = length(Y00{n});
    Y(n) = mean(Y00{n});
    stdY(n) = std(Y00{n});
end


switch xIncremenScale
    case 'log'
        X = exp(X);
end

end
