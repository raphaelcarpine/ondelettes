function plotSpectrums(Freqs, Spectrums, spectrumFrequencyScale, spectrumScale, figName)
%PLOTSPECTRUMS Summary of this function goes here
%   Detailed explanation goes here

fig = figure('Name', figName);
ax = axes(fig);
plts = plot(ax, Freqs, Spectrums);
set(ax, 'XScale', spectrumFrequencyScale, 'YScale', spectrumScale);
xlabel(ax, 'Frequency [Hz]');
ylabel(ax, 'Amplitude CWT');

% legend
shockTimeNames = cell(size(Spectrums, 1), 1);
for k_shock = 1:size(Spectrums, 1)
    shockTimeNames{k_shock} = ['t_{', num2str(k_shock), '}'];
end

if length(shockTimeNames) <= 10
    legend(ax, shockTimeNames, 'Location', 'best');
end

% data tips
for k_plot = 1:length(plts)
    plt = plts(k_plot);
    plt.DataTipTemplate.DataTipRows = [plt.DataTipTemplate.DataTipRows(1); plt.DataTipTemplate.DataTipRows];
    plt.DataTipTemplate.DataTipRows(1).Label = shockTimeNames{k_plot};
    plt.DataTipTemplate.DataTipRows(1).Value = '';
    plt.DataTipTemplate.DataTipRows(2).Label = 'freq';
    plt.DataTipTemplate.DataTipRows(3).Label = 'ampl';
end

end

