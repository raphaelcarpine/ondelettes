function map = CWT(t, X, freqs, Nt, ondelette, Q)
%CWT Summary of this function goes here
%   Detailed explanation goes here

% p = inputParser ;
% % parametres par defaut
% ZeroPaddingDef = 1;
% ctDef = 5;
% addRequired(p,'X')
% addRequired(p,'Y')
% addRequired(p,'WvltFreq')
% addRequired(p,'Q')
% addParameter(p,'ZeroPadding',ZeroPaddingDef);
% addParameter(p,'ct',ctDef);
% 
% parse(p,X,WvltFreq,Y,Qin,varargin{:});
% 
% %
% ZeroPadding = p.Results.ZeroPadding;
% ct = p.Results.ct;

switch nargin
    case 4
        ondelette = 'gabor';
        Q = 1;
    case 5
        Q = 1;
end



ondelettes = containers.Map;
ondelettes('gabor') = @(Q) @(t) exp(2*1i*pi*t - t.^2/(Q^2/pi^2));
ondelettes('cauchy') = @(Q) @(t) (1i ./ (2*pi*t/(2*Q^2) + 1i)).^(2*Q^2+1/2);
ondelettes('fourier') = @(Q) @(t) 1/sqrt(2*pi)*exp(2*1i*pi*t);


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
    freqs2(:, j) = freqs';
end
for i = 1:length(freqs)
    T2(i, :) = T;
end

figure;
surf(freqs2, T2, log(abs(map)), 'edgecolor', 'none');
xlabel('f');
ylabel('t');
zlabel('log(arg)');


% phi = angle(map);
% for ifreq = 1:length(freqs)
%     for iT = 1:length(T)-1
%         delta = phi(ifreq, iT+1) - phi(ifreq, iT);
%         if abs(delta + 2*pi) < abs(delta)
%             phi(ifreq, iT+1:end) = phi(ifreq, iT+1:end) + 2*pi;
%         elseif abs(delta - 2*pi) < abs(delta)
%             phi(ifreq, iT+1:end) = phi(ifreq, iT+1:end) - 2*pi;
%         end
%     end
% end
% iT = int16(length(T)/2);
% for ifreq = 1:length(freqs)-1
%     delta = phi(ifreq+1, iT) - phi(ifreq, iT);
%     if abs(delta + 2*pi) < abs(delta)
%         phi(ifreq+1:end, :) = phi(ifreq+1:end, :) + 2*pi;
%     elseif abs(delta - 2*pi) < abs(delta)
%         phi(ifreq+1:end, :) = phi(ifreq+1:end, :) - 2*pi;
%     end
% end
% figure;
% surf(freqs2, T2, phi + 2*pi*freqs'*T, 'edgecolor', 'none');
% xlabel('f');
% ylabel('t');
% zlabel('phi+omega*t');

end