function map = WvltCompTemp(t, X, freqs, Nt, Q)
%CWT Summary of this function goes here
%   Detailed explanation goes here


n = 2*Q^2-1/2;
psi = @(t) exp(log(n+1) * 1i ./ (2*pi*t/(n+1/2) + 1i));

T = linspace(t(1), t(end), Nt);
  
map = zeros(length(freqs), length(T));
wait = waitbar(0, sprintf('0/%d', length(freqs)), 'Name', "Calcul de la transformée");
for ifreq = 1:length(freqs)
    for iT = 1:length(T)
        map(ifreq, iT) = trapz(t, X.*conj(psi(freqs(ifreq)*(t-T(iT)))));
    end
    waitbar(ifreq/length(freqs), wait, sprintf('%d/%d',ifreq,length(freqs)));
end
close(wait);

end