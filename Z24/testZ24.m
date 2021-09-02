clear all

measure = 1; % 1-17
config = 5; % 1-9

[X, t, labels] = getDataZ24(measure, config);

%% 

fig = figure('Name', sprintf('measure %02u, config %u', [measure, config]));
plt = plot(t, X);
xlabel('Time [s]');
ylabel('Acceleratio [m/s²]');
legend(labels, 'NumColumns', 2);
selectLine(gca);

% CWT

fmin = 0.1;
fmax = 10;
Q = 10;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q);

%% mode 1

I = false(size(labels));
for ki = 1:length(I)
    if labels{ki}(1) ~= 'R' && labels{ki}(end) == 'V'
        I(ki) = true;
    end
end

x = sum(X(I, :), 1);


fig = figure('Name', sprintf('measure %02u, config %u', [measure, config]));
plt2 = plot(t, x);
xlabel('Time [s]');
ylabel('Acceleratio [m/s²]');
legend(labels, 'NumColumns', 2);
selectLine(gca);

% CWT

fmin = 0.1;
fmax = 10;
Q = 10;

WaveletMenu('WaveletPlot', plt2, 'fmin', fmin, 'fmax', fmax, 'Q', Q);