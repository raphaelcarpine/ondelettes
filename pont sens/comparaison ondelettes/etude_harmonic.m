% etude des effets de bord des ondelettes harminuqe et littlewood-paley
eulerGamma = 0.5772156649;

%%
N = 7;

x = linspace(-N*pi, N*pi, 1001);
sinintX = nan(size(x));
cosintX = nan(size(x));
for ix = 1:length(x)
    sinintX(ix) = sinint(x(ix));
    cosintX(ix) = cosint(x(ix));
end

figure;
hold on
plot(x, cosintX);
plot(x, sinintX);
legend({'Ci(x) = Re\{Eii(x)\}', 'Si(x) = Im\{Eii(x)\}'}, 'Location', 'Best');
xlim([x(1), x(end)]);
ylim([-2 2]);
Xt = [0];
Xtl = {'0'};
n = 2;
while n <= N
    Xt = [-n*pi, Xt, n*pi];
    Xtl = cat(2, {['-', num2str(n), '\pi']}, Xtl);
    Xtl = cat(2, Xtl, {[num2str(n), '\pi']});
    n = n+2;
end
xticks(Xt);
xticklabels(Xtl);
yticks([-pi/2, 0, pi/2]);
yticklabels({'-\pi/2', '0', '\pi/2'});

%%

% Lambda = @(x) 1/(2*pi) * (sign(x)*(pi/2-sinint(abs(x)))  + 1i*cosint(abs(x)));
% 
% x = linspace(-5*pi, 5*pi, 1001);
% Lx = nan(size(x));
% for ix = 1:length(x)
%     Lx(ix) = Lambda(x(ix));
% end
% figure;
% plot(x, real(Lx));
% hold on
% plot(x, imag(Lx));

%% Lambda(a, b) en fonction de 1/a

alpha = 1;
beta = 2;
b = 5;
nu0 = 4;
F = linspace(0, 5, 1000); % 1/a
F = [F, linspace(1.9, 2.1, 1000)];
F = [F, linspace(3.9, 4.1, 1000)];
F = sort(F);

R = @(f) sign(beta*f-nu0)*pi/2 + 1i*Eii(2*pi*(beta*f-nu0)*b) - ( sign(alpha*f-nu0)*pi/2 + 1i*Eii(2*pi*(alpha*f-nu0)*b) );

Rf = nan(size(F));
for kf = 1:length(F)
    Rf(kf) = R(F(kf));
end

rRfnan = real(Rf);
Fnan = F;
kf = 1;
while kf < length(Fnan)
    if abs(rRfnan(kf+1) - rRfnan(kf)) > pi/2
        rRfnan = [rRfnan(1:kf), nan, rRfnan(kf+1:end)];
        Fnan = [Fnan(1:kf), (Fnan(kf)+Fnan(kf+1))/2, Fnan(kf+1:end)];
        kf = kf+2;
    else
        kf = kf+1;
    end
end

fig = figure;
fig.Position = fig.Position .* [1 1 0 0] + [0 0 420 300];
hold on
Rline = plot(Fnan, rRfnan);
Iline = plot(F, imag(Rf));
Absline = plot(F, abs(Rf), 'Color', 0.7*[1 0 0]);
ax = gca;
ax.ColorOrderIndex = 1;
Rline2 = plot(F, real(Rf), '--');
% uistack(Rline2, 'top');
% uistack(Rline, 'top');
% uistack(Absline, 'top');
xlim([F(1), F(end)]);
ylim([-3 3]);
xlabel('1/a');
legend({'Re\{\Lambda(a, b)\}', 'Im\{\Lambda(a, b)\}', '|\Lambda(a, b)|'}, 'Location', 'NorthWest');
xticks([0, nu0/beta, nu0/alpha]);
xticklabels({'0', '\nu_0 /\beta', '\nu_0 /\alpha'});
yticks([-pi/2, 0, pi/2]);
yticklabels({'-\pi/2', '0', '\pi/2'});

%% Lambda(a, b) en fonction de b

alpha = 1;
beta = 2;
a = 1/1.5;
nu0 = 4;
TauAlpha = 1/abs(alpha/a-nu0);
B = linspace(-10*TauAlpha, 10*TauAlpha, 1000);

R = @(b) sign(beta/a-nu0)*pi/2 + 1i*Eii(2*pi*(beta/a-nu0)*b) - ( sign(alpha/a-nu0)*pi/2 + 1i*Eii(2*pi*(alpha/a-nu0)*b) );

Rb = nan(size(B));
for kb = 1:length(B)
    Rb(kb) = R(B(kb));
end

fig = figure;
fig.Position = fig.Position .* [1 1 0 0] + [0 0 420 300];
hold on
Rline = plot(B, real(Rb));
Iline = plot(B, imag(Rb));
Absline = plot(B, abs(Rb), 'Color', 0.7*[1 0 0]);
% uistack(Rline2, 'top');
% uistack(Rline, 'top');
% uistack(Absline, 'top');
xlim([B(1), B(end)]);
xlabel('b');
legend({'Re\{\Lambda(a, b)\}', 'Im\{\Lambda(a, b)\}', '|\Lambda(a, b)|'}, 'Location', 'NorthWest');
Xt = [0];
Xtl = {'0'};
n = 2;
while n <= N
    Xt = [-n*pi, Xt, n*pi];
    Xtl = cat(2, {['-', num2str(n), '\pi']}, Xtl);
    Xtl = cat(2, Xtl, {[num2str(n), '\pi']});
    n = n+2;
end
% xticks([0, nu0/beta, nu0/alpha]);
% xticklabels({'0', '\nu_0 /\beta', '\nu_0 /\alpha'});
% yticks([-pi/2, 0, pi/2]);
% yticklabels({'-\pi/2', '0', '\pi/2'});

%% Lambda(a, b) en fonction de 1/a et b

alpha = 1;
beta = 2;
nu0 = 4;
TauAlpha = 1/abs(alpha/a-nu0);

F = linspace(0, 6, 100);
F = [F, 2 + logspace(-8, 0, 30)];
F = [F, 2 - logspace(-8, 0, 30)];
F = [F, 4 + logspace(-8, 0, 30)];
F = [F, 4 - logspace(-8, 0, 30)];
F = sort(F);

B = linspace(-10, 10, 100);
B = [B, 0 + logspace(-8, -1, 30)];
B = [B, 0 - logspace(-8, -1, 30)];
B = sort(B);

B = linspace(0, 10, 100);
B = [B, 0 + logspace(-8, -1, 50)];
B = sort(B);

R = @(f, b) sign(beta*f-nu0)*pi/2 + 1i*Eii(2*pi*(beta*f-nu0)*b) - ( sign(alpha*f-nu0)*pi/2 + 1i*Eii(2*pi*(alpha*f-nu0)*b) );

Rfb = nan(length(F), length(B));
for kf = 1:length(F)
    for kb = 1:length(B)
        Rfb(kf, kb) = R(F(kf), B(kb));
    end
end

[Bgrid, Fgrid] = meshgrid(B, F);

%%

figure;
surf(Fgrid, Bgrid, abs(Rfb), 'EdgeColor', 'none', 'FaceColor', [0 0.4470 0.7410], 'FaceLighting', 'gouraud');
xlabel('1/a');
ylabel('b');
zlabel('|\Lambda(a,b)|'),
zlim([0, 8]);
% zticks([0, 2*pi]);
% zticklabels({'0', '2\pi'});
xticks([0, nu0/beta, nu0/alpha]);
xticklabels({'0', '\nu_0/\beta', '\nu_0/\alpha'});
set(gca, 'XDir','reverse')
light();

%% T[U](a, b)

alpha = 1;
beta = 2;
nu0 = 4;

F = linspace(0, 6, 100);
F = [F, 2 + logspace(-8, 0, 30)];
F = [F, 2 - logspace(-8, 0, 30)];
F = [F, 4 + logspace(-8, 0, 30)];
F = [F, 4 - logspace(-8, 0, 30)];
F = sort(F);

B = linspace(0, 10, 2);

TUfb = nan(length(F), length(B));
for kf = 1:length(F)
    for kb = 1:length(B)
        TUfb(kf, kb) = (F(kf) >= nu0/beta) & (F(kf) <= nu0/alpha);
    end
end

[Bgrid, Fgrid] = meshgrid(B, F);

figure;
surf(Fgrid, Bgrid, abs(TUfb), 'EdgeColor', 'none', 'FaceColor', [0 0.4470 0.7410], 'FaceLighting', 'gouraud');
xlabel('1/a');
ylabel('b');
zlabel('|T[U](a,b)|'),
zlim([0, 1.5]);
% zticks([0, 2*pi]);
% zticklabels({'0', '2\pi'});
xticks([0, nu0/beta, nu0/alpha]);
xticklabels({'0', '\nu_0/\beta', '\nu_0/\alpha'});
zticks([0, 1]);
set(gca, 'XDir','reverse')
light();


%% fonctions

function I = sinint(x)
if x == 0
    I = 0;
    return
end

dx = 0.01;

X = linspace(0, x, ceil(abs(x/dx)) + 1);
dx = mean(diff(X));

Y = sinc(X/pi);

I = (sum(Y) - (Y(1)+Y(end))/2) * dx;
end

function I = cosint(x)
x = abs(x);

if x == 0
    I = nan;
    return
end

dx = 0.01;

X = linspace(0, x, ceil(x/dx) + 1);
dx = mean(diff(X));

Y = (cos(X) - 1) ./ X;
Y(1) = 0;

I = (sum(Y) - (Y(1)+Y(end))/2) * dx;

I = 0.5772156649 + log(x) + I;
end

function I = Eii(x)
I = cosint(x) + 1i*sinint(x);
end