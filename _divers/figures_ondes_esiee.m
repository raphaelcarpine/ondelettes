f = 1;
phi = 1;
L = 1;
zeta = 0.03;
x = linspace(0, 10, 1000);
Nframes = 20;

%%

gifName = 'onde_amortie.gif';


%%
k = 2*pi/L * (sqrt(1-zeta^2) - 1i*zeta);
w = 2*pi*f;
U = @(t) real(exp(1i*((w*t - k*x + phi))));

t = 1/f * (0:Nframes-1)/Nframes;

figure();
plt = plot(x, U(t(1)));
ylim([-1.1, 1.1]);
xlabel('x');
ylabel('U');
xticks([]);
yticks([]);

set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');

input(' ');


gif(gifName, 'frame', gcf, 'DelayTime', mean(diff(t)));

for kt = 2:Nframes
    plt.YData = U(t(kt));
    drawnow;
    gif;
end







