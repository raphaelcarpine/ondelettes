clear all

file1 = 'maquetteTMDesiee2020\data\20201027_TMD_converted\mData.mat';
file2 = 'maquetteTMDesiee2020\data\20201112_TMD_converted\mData.mat';

%% TMD seul amorti

load(file1);

k_essai = 2; % [1, ..., 10]

tX = TMDseulAmort;
tX  = transpose(tX);
t = tX(1, :);
X = tX(2, :);

ti = [4.4, 37.3, 57.1, 73, 85.4, 96, 105, 113.5, 119.8, 126];
tf = [36.8, 56.7, 71.6, 85, 96, 105, 113.5, 119.8, 125];


% plot
fig = figure;
ax = axes(fig);
plt = plot(ax, t, X);
xlabel(ax, 'Time [s]');
ylabel(ax, 'Acceleration [m/s²]');


% CWT
fmin = 0.3;
fmax = 1.5;
Q = 4;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
    'MaxParallelRidges', 1, 'MaxRidges', 1, 'XLim', [ti(k_essai), tf(k_essai)],...
    'WvltScale', 'lin', 'RemoveMean', true);

% "time, ampl  |lin|  |log|" -> "time, ampl  |lin|  |lin|"


%% TMD seul non amorti

load(file2);

tX = TMDseulnonamorti_v2;
tX  = transpose(tX);
t = tX(1, :);
X = tX(2, :); % acceleration ch1

ti = 3;
tf = 26;

% plot
fig = figure;
ax = axes(fig);
plt = plot(ax, t, X);
xlabel(ax, 'Time [s]');
ylabel(ax, 'Acceleration [m/s²]');


% CWT
fmin = 0.3;
fmax = 1.5;
Q = 10;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
    'MaxParallelRidges', 1, 'MaxRidges', 1, 'XLim', [ti, tf],...
    'WvltScale', 'lin', 'RemoveMean', true);

% "time, ampl  |lin|  |log|" -> "time, ampl  |lin|  |lin|"



%% maquette + TMD non amorti

load(file2);
tX = MaquetteplusTMDnonamortie_V2;
ti = 0;
tf = 30;

% load(file1);
% tX = TMDplusMaquetteNonAmort_2978ch2;
% ti = 7.74;
% tf = 92;

tX  = transpose(tX);
t = tX(1, :);
X = tX(2, :); % acceleration ch1


% plot
fig = figure;
ax = axes(fig);
plt = plot(ax, t, X);
xlabel(ax, 'Time [s]');
ylabel(ax, 'Acceleration [m/s²]');

% CWT
fmin = 0.3;
fmax = 1.5;
Q = 5;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
    'MaxParallelRidges', 1, 'MaxRidges', 1, 'XLim', [ti, tf],...
    'WvltScale', 'lin', 'RemoveMean', true);

% "time, ampl  |lin|  |log|" -> "time, ampl  |lin|  |lin|"



%% maquette + TMD amorti

load(file2);
tX = MaquetteplusTMDamortie_V3;
ti = 15.5;
tf = 96;
k_essai = 1;

% load(file1);
% tX = TMDplusMaquetteAmort_2977ch3;
% ti = [9.73, 51.1, 77.8];
% tf = [40, 74, 102];
% k_essai = 1;

tX  = transpose(tX);
t = tX(1, :);
X = tX(2, :); % acceleration ch1

% plot
fig = figure;
ax = axes(fig);
plt = plot(ax, t, X);
xlabel(ax, 'Time [s]');
ylabel(ax, 'Acceleration [m/s²]');

% CWT
fmin = 0.3;
fmax = 1.5;
Q = 5;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
    'MaxParallelRidges', 1, 'MaxRidges', 1, 'XLim', [ti(k_essai), tf(k_essai)],...
    'WvltScale', 'lin', 'RemoveMean', true);

% "time, ampl  |lin|  |log|" -> "time, ampl  |lin|  |lin|"