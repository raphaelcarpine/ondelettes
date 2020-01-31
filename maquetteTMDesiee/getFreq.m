function f = getFreq(t, x, fbounds)
%GETFREQ Summary of this function goes here
%   Detailed explanation goes here

dt = (t(end)-t(1))/(length(t)-1);

if max(abs(diff(t)/dt-1)) > 1e-5
    warning(['pas de temps non constant, erreur : ', num2str(max(abs(diff(t)/dt-1)))]);
end


% fft
T = t(end)-t(1);
fx = abs(fft(x));
freqs = 1:length(fx);
freqs = (freqs-1) / T;

% bounds
fx = fx(freqs>=fbounds(1) & freqs<=fbounds(2));
freqs = freqs(freqs>=fbounds(1) & freqs<=fbounds(2));

% maximum
[~, kfmax] = max(fx);
f = freqs(kfmax);

% fig = figure;
% ax = axes(fig);
% hold(ax, 'on');
% plot(ax, freqs, fx);
% plot(ax, f, max(fx), 'r+');

end

