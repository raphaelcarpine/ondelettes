d = 6;
L = 5;
N = 10;
c = 1;


k = 3;

phi_k = @(x) sin(k*pi*x/L);

sndMb = @(t) sum( phi_k((0:N-1)*d + c*t) .* (0 <= ((0:N-1)*d + c*t) & ((0:N-1)*d + c*t) <= L) );


T = linspace(-(L+N*d)/c, 2*L/c, 10000);
SndMb = nan(size(T));
for it = 1:length(T)
    SndMb(it) = sndMb(T(it));
end


figure;
plot(T, SndMb);
xlabel('t');
xticks([]);
ylabel('$\langle\varphi_k, \mu_q(t)\rangle$', 'Interpreter', 'latex');
yticks([]);




