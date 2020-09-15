beta = 1;
delta = 0.5;
nu0 = 4;

% maxR = @(a, b) 1/2*(1-erf(b/(a*delta*sqrt(2))))*exp(-2*pi^2*delta^2*(beta-a*nu0)^2) + 1/2*erf(pi*delta*sqrt(2)*(beta-a*nu0))*exp(-b^2/(2*a^2*delta^2));
% maxR = @(a, b) 1/2*(1-erf(b/(a*delta*sqrt(2))));


R(1/4, 0.1, beta, delta, nu0);


F = linspace(0, 12, 100);
F = [F, 4 + logspace(-8, 0, 30)];
F = [F, 4 - logspace(-8, 0, 30)];
F = sort(F);

B = linspace(-0.5, 0.5, 100);
B = [B, 0 + logspace(-8, -1, 30)];
B = [B, 0 - logspace(-8, -1, 30)];
B = sort(B);
B = B(B>=0);

% B = linspace(0, 0.5, 100);
% B = [B, 0 + logspace(-8, -1, 50)];
% B = sort(B);

Rf = @(f, b) R(1/f, b, beta, delta, nu0);

Rfb = nan(length(F), length(B));
for kf = 1:length(F)
    for kb = 1:length(B)
        Rfb(kf, kb) = Rf(F(kf), B(kb));
    end
end

[Bgrid, Fgrid] = meshgrid(B, F);

%%

figure;
surf(Fgrid, Bgrid, abs(Rfb), 'EdgeColor', 'none', 'FaceColor', [0 0.4470 0.7410], 'FaceLighting', 'gouraud');
xlabel('1/a');
ylabel('b');
zlabel('|R(a,b)|'),
zlim([0, 0.6]);
xlim([F(1), F(end)]);
ylim([B(1), B(end)]);
% zticks([0, 2*pi]);
% zticklabels({'0', '2\pi'});
xticks([0, nu0/beta]);
xticklabels({'0', '\nu_0/\beta'});
set(gca, 'YDir','reverse')
light();

%% T

B = [0, 10];

Tu = nan(length(F), length(B));
for kf = 1:length(F)
    for kb = 1:length(B)
        Tu(kf, kb) = exp(-2*pi^2 * delta^2 * (nu0./F(kf) - beta).^2);
    end
end

[Bgrid, Fgrid] = meshgrid(B, F);

figure;
surf(Fgrid, Bgrid, Tu, 'EdgeColor', 'none', 'FaceColor', [0 0.4470 0.7410], 'FaceLighting', 'gouraud');
xlabel('1/a');
ylabel('b');
zlabel('|T[U](a,b)|'),
zlim([0, 1.3]);
xlim([F(1), F(end)]);
ylim([B(1), B(end)]);
zticks([0, 1]);
xticks([0, nu0/beta]);
xticklabels({'0', '\nu_0/\beta'});
set(gca, 'YDir','reverse')
light();












%%


function I = R(a, b, beta, delta, nu0)
N = 5;
n = 200;

dt = 1 / (n*abs(beta-a*nu0));
dt = min(dt, N*delta/1000);

t1 = max(-N*delta, b/a);
t2 = b/a*(b/a>=0) + N*delta;
t = t1:dt:t2;

I = 1/(delta*sqrt(2*pi)) * sum(exp(-t.^2/(2*delta^2) + 2i*pi*(beta-a*nu0)*t)) * dt;

end