%% loading

folder = 'simulations michalis\data';
n = 4;


X = nan(0);
stories = {'first', 'second'};
for story = 1:2
    fileName = ['abs_acc_', stories{story}, 'StoryAthens_', num2str(n)];
    S = load(fullfile(folder, fileName));
    field = fieldnames(S);
    
    X = [X; transpose(S.(field{1}))];
end


%% plot

Fs = 100; % sampling frequency
T = (0:size(X, 2)-1) / Fs;


figure;
plt = plot(T, X);
xlabel('Time [t]');
ylabel('Acceleration [m/s²]');
legend({'1st story', '2nd story'});


%% wavelet

fmin = 0.5;
fmax = 10;
% fmin = 10;
% fmax = 50;
Q = 5;
Multi = true;

WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q,...
    'MultiSIgnalMode', Multi);




