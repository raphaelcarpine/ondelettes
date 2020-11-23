function plotSpectrums(Freqs, Spectrums, spectrumFrequencyScale, spectrumScale, figName)
%PLOTSPECTRUMS Summary of this function goes here
%   Detailed explanation goes here

fig = figure('Name', figName);
ax = axes(fig);
plot(ax, Freqs, Spectrums);
set(ax, 'XScale', spectrumFrequencyScale, 'YScale', spectrumScale);
xlabel(ax, 'Frequency [Hz]');
ylabel(ax, 'Amplitude CWT');

shockTimeNames = ['t_{'.*ones(size(Spectrums, 1), 1), num2str(transpose(1:size(Spectrums, 1))), '}'.*ones(size(Spectrums, 1), 1)];
shockTimeNames = cellstr(shockTimeNames);

legend(ax, shockTimeNames, 'Location', 'best');
end

