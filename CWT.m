function CWT(t, X, freqs, Nt, ondelette, Q)
%CWT Summary of this function goes here
%   Detailed explanation goes here
switch nargin
    case 4
        ondelette = 'gabor';
        Q = 1;
    case 5
        Q = 1;
end

ondelettes = containers.Map;
ondelettes('gabor') = @(Q) @(t) exp(2*1i*pi*t - t.^2*Q);


Psi = ondelettes(ondelette);
psi = Psi(Q);

T = linspace(t(1), t(end), Nt);
  
map = zeros(length(freqs), length(T));
wait = waitbar(0, sprintf('0/%d', length(freqs)), 'Name', "Calcul de la transformée");
for ifreq = 1:length(freqs)
    for iT = 1:length(T)
        map(ifreq, iT) = trapz(t, X.*psi(freqs(ifreq)*(t-T(iT))));
    end
    waitbar(ifreq/length(freqs), wait, sprintf('%d/%d',ifreq,length(freqs)));
end
close(wait);

for j = 1:length(T)
    freqs2(:, j) = freqs;
end
for i = 1:length(freqs)
    T2(i, :) = T.';
end

figure;
surf(freqs2, T2, abs(map), 'edgecolor', 'none');
xlabel('f');
ylabel('t');

end