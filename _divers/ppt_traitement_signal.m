f0 = 10;
zeta = 0.01;
A = 1;
T = [-1, 100];
Fe = 500;

lambd = 2*pi*f0*(1i*sqrt(1-zeta^2) - zeta);
t = T(1):1/Fe:T(2);
s = imag(A*(t>=0).*exp(lambd*t));

fig = figure;
plot(t, s);
xlabel('t');
ylabel('s(t)');
xlim([-1 10]);
fig.Position(3:4) = [480 320];

S = fft(s);
S = S/Fe;
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

S = S.*exp(2i*pi*f); % translation

fig = figure;
plot(f, abs(S));
xlabel('\nu');
ylabel('|S(\nu)|');
xlim([-20 20]);
ylim([0 1]);
fig.Position(3:4) = [480 320];

fig = figure;
plot(f, angle(S));
xlabel('\nu');
ylabel('arg S(\nu)');
xlim([-20 20]);
% ylim([-pi pi]);
fig.Position(3:4) = [480 320];



%%

f0 = pi;
A = 1;
T = [-5, 1000];
Fe = 500;

lambd = 2i*pi*f0
t = T(1):1/Fe:T(2);
s = imag(A*exp(lambd*t));
s2 = s;
s2(~(t >= 0 & t <= 3)) = nan;

fig = figure;
plot(t, s, '--');
hold on
set(gca,'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
plot(t, s2);
xlabel('t');
ylabel('s(t)');
xlim([-1.5 6]);
ylim(1.2*[-1 1]);
fig.Position(3:4) = [480 320];
legend({'s(t)', 'š(t)'});

S = fft(s);
S = S/Fe;
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

S = S.*exp(2i*pi*f); % translation

fig = figure;
plot(f, abs(S), '--');
xlabel('\nu');
ylabel('|S(\nu)|');
xlim(10*[-1 1]);
% ylim([0 1]);
yticks([]);
fig.Position(3:4) = [480 320];



%%

f0 = pi;
A = 1;
T = [-5, 1000];
Fe = 500;

lambd = 2i*pi*f0;
t = T(1):1/Fe:T(2);
s = imag(A*exp(lambd*mod(t, 3)));

fig = figure;
plot(t, s);
xlabel('t');
ylabel('s(t) (périodisé)');
xlim([-1.5 10]);
ylim(1.2*[-1 1]);
fig.Position(3:4) = [480 320];

S = fft(s);
S = S/Fe;
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

S = S.*exp(2i*pi*f); % translation

fig = figure;
plot(f, abs(S));
xlabel('\nu_k');
ylabel('|S_k|');
xlim(10*[-1 1]);
% ylim([0 1]);
fig.Position(3:4) = [480 320];
yticks([]);



%%

f0 = pi;
A = 1;
T = [-5, 1000];
Fe = 500;

lambd = 2i*pi*f0;
t = T(1):1/Fe:T(2);
s = (t>=0 & t<=3).*imag(A*exp(lambd*t));

fig = figure;
plot(t, s);
xlabel('t');
ylabel('s(t) (fenêtré))');
xlim([-1.5 10]);
ylim(1.2*[-1 1]);
fig.Position(3:4) = [480 320];

S = fft(s);
S = S/Fe;
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

fig = figure;
plot(f, abs(S));
xlabel('\nu');
ylabel('|S(\nu)|');
xlim(10*[-1 1]);
% ylim([0 1]);
fig.Position(3:4) = [480 320];
yticks([]);


t = 0:1/Fe:3;
s = imag(A*exp(lambd*t));

S = fft(s);
S = S/Fe;
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

hold on
set(gca,'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter(f, abs(S), 'o');

%%

Fe = 500;
t = -100:1/Fe:100;

s1 = t>=0 & t<1; % porte
s2 = s1.*(0.5-0.5*cos(2*pi*t)); % hanning
s3 = s1.*(0.54-0.46*cos(2*pi*t)); % hamming
s4 = s1.*(0.42-0.5*cos(2*pi*t)+0.08*cos(4*pi*t)); % blackman

% s2 = s2/0.75;
% s3 = s3/(0.56+0.23);
% s4 = s4/(0.42+0.25+0.04);

fig = figure;
plot(t, s1);
hold on
plot(t, s2);
plot(t, s3);
plot(t, s4);
xlabel('t');
ylabel('u(t)');
xlim([-0.3 1.3]);
ylim([-0.1 1.1]);
xticks([]);
yticks([0 1]);
fig.Position(3:4) = [480 320];
legend({'porte', 'Hanning', 'Hamming', 'Blackman'});





S = fft(s1);
S = S/Fe;
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

fig = figure;
plot(f, abs(S));
hold on
xlabel('\nu');
ylabel('|U(\nu)|');
xlim(10*[-1 1]);
ylim([0.00001 1.5]);
fig.Position(3:4) = [480 320];
set(gca, 'YScale', 'log');
xticks([]);



S = fft(s2);
S = S/Fe;
S = S/S(1);
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

plot(f, abs(S));



S = fft(s3);
S = S/Fe;
S = S/S(1);
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;
plot(f, abs(S));



S = fft(s4);
S = S/Fe;
S = S/S(1);
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

plot(f, abs(S));

legend({'porte', 'Hanning', 'Hamming', 'Blackman'});

%%

f0 = 10;
A = 1;
T = [0, 10000];
Fe = 500;
Fe2 = 154;

lambd = 2i*pi*f0;
t = T(1):1/Fe:T(2);
s = imag(A*exp(lambd*t));
t2 = T(1):1/Fe2:T(2);
s2 = imag(A*exp(lambd*t2));

fig = figure;
plot(t, s);
hold on
set(gca,'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter(t2, s2, 'r');
xlabel('t');
ylabel('s');
xlim([0 1]);
ylim(1.1*[-1 1]);
fig.Position(3:4) = [480 320];
legend({'s(t)', 's_k'});

%%

f0 = 10;
A = 1;
T = [0, 10000+pi];
Fe2 = 17;
Fe = 400;

lambd = 2i*pi*f0;
t = T(1):1/Fe:T(2);
s = imag(A*exp(lambd*t));
t2 = T(1):1/Fe2:T(2);
s2 = imag(A*exp(lambd*t2));

fig = figure;
plot(t, s);
hold on
set(gca,'ColorOrderIndex', get(gca,'ColorOrderIndex')-1);
scatter(t2, s2, 'r');
xlabel('t');
ylabel('s');
xlim([0 1]);
ylim(1.1*[-1 1]);
fig.Position(3:4) = [480 320];
legend({'s(t)', 's_k'});
xticks([]);
yticks([]);

s3 = -sin(2*pi*(Fe2-f0)*t);

plot(t, s3, 'r--');


S = fft(s);
S = S/Fe;
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
I = max(abs(S));
I = find(abs(S) == I);
S(I) = 5000;
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

fig = figure;
plot(f, abs(S), 'LineWidth', 2);
hold on
xlabel('\nu');
ylabel('|S|');
xlim(15*[-1 1]);
ylim([0 6000]);
fig.Position(3:4) = [480 320];
xticks([-Fe2/2, 0, Fe2/2]);
xticklabels({'-F_e/2', '0', 'F_e/2'});
yticks([]);



S = fft(s2);
S = S/Fe2;
N = length(S);
S = S([floor(N/2)+1:end, 1:floor(N/2)]);
I = max(abs(S));
I = find(abs(S) == I);
S(I) = 5000;
Df = Fe2/length(t2);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

plot(f, abs(S), 'r', 'LineWidth', 2);

%%

t = 0:0.001:10;
s = zeros(size(t));
s(1) = 1;

figure;
plt = plot(t, s);

WaveletMenu('WaveletPlot', plt);

%%

f0 = 10;
zeta = 0.01;
Fe = 400;

syst = systLin(1, (2*pi*f0)^2, 2*zeta*2*pi*f0);

dt = 1/Fe
t = 0:dt:40;
f = randn(size(t));

s = syst.response(f, dt, 1);

fig = figure;
plot(t, s);
xlabel('t');
ylabel('s');
% xlim([0 1]);
ylim(max(abs(s))*1.1*[-1 1]);
fig.Position(3:4) = [350 220];
% legend({'s(t)', 's_k'});
xticks([]);
yticks([]);

%%

f0 = 10;
T = [0, 0.5];
Fe = 400;

lambd = 2i*pi*f0;
t = T(1):1/Fe:T(2);
sigma = 0.1;

s1 = imag(exp(lambd*t)) + sigma*randn(size(s));
s2 = imag(exp(lambd*t + 0.65i*pi)) + sigma*randn(size(s));

d = 2.5;
figure;
plot(t, s1);
hold on
plot(t, s2+d);
yline(0);
yline(d);
ylim([-1.5, d+1.5]);
axis off


%%

f0 = 10;
zeta = 0.03;
A = 1;
T = [0, 100];
Fe = 500;
sigma = 0.05;

for k = 1:4
    
    lambd = 2*pi*f0*(1i*sqrt(1-zeta^2) - zeta);
    t = T(1):1/Fe:T(2);
    s = imag(A*(t>=0).*exp(lambd*t));
    s = s + sigma*randn(size(s));
    
    
    fig = figure;
    plot(t, s);
    xlabel('t');
    ylabel('s(t)');
    xlim([0 4]);
    ylim(1.1*[-1 1]);
    yticks([]);
    fig.Position(3:4) = [350 220];
    
    S = fft(s);
    S = S/Fe;
    N = length(S);
    S = S([floor(N/2)+1:end, 1:floor(N/2)]);
    Df = Fe/length(t);
    f = (-N+floor(N/2):floor(N/2)-1) * Df;
    
    fig = figure;
    yyaxis left
    plot(f, abs(S));
    xlabel('\nu');
    ylabel('|S(\nu)|');
    yticks(0);
    yyaxis right
    plot(f, angle(S*exp(0.5i*pi))-0.5*pi);
    ylabel('arg S(\nu)');
    yticks([-pi, 0]);
    yticklabels({'-\pi', '0'});
    % ylim([-pi 0] + 0.*[-1 1]);
    xlim([0 17]);
    % ylim([0 1]);
    fig.Position(3:4) = [350 220];
end



s0 = nan(0, length(s));
S0 = nan(0, length(S));

for k = 1:200
    s = imag(A*(t>=0).*exp(lambd*t));
    s = s + sigma*randn(size(s));
    
    s0(end+1, :) = s;
    
    S = fft(s);
    S = S/Fe;
    N = length(S);
    S = S([floor(N/2)+1:end, 1:floor(N/2)]);
    
    S0(end+1, :) = S;
end

s = mean(s0, 1);
S = mean(S0, 1);


fig = figure;
plot(t, s);
xlabel('t');
ylabel('s(t)');
xlim([0 4]);
ylim(1.1*[-1 1]);
yticks([]);
fig.Position(3:4) = [480 320];

fig = figure;
yyaxis left
plot(f, abs(S));
xlabel('\nu');
ylabel('|S(\nu)|');
yticks([0]);
yyaxis right
plot(f, angle(S*exp(0.5i*pi))-0.5*pi);
ylabel('arg S(\nu)');
yticks([-pi, 0]);
yticklabels({'-\pi', '0'});
% ylim([-pi 0] + 0.*[-1 1]);
xlim([0 17]);
% ylim([0 1]);
fig.Position(3:4) = [480 320];


%%

f0 = 10;
zeta = 0.03;
T = [0, 100];
Fe = 500;
sigma = 0.05;
F0 = 1000;

for k = 1:3
    syst = systLin(1, (2*pi*f0)^2, 2*zeta*2*pi*f0);
    
    dt = 1/Fe;
    t = 0:dt:100;
    f = F0*randn(size(t));
    
    s = syst.response(f, dt, 1);
    
    
    fig = figure;
    plot(t, s);
    xlabel('t');
    ylabel('s(t)');
    xlim([0 4]);
    ylim(1.1*[-1 1]);
    yticks([]);
    fig.Position(3:4) = [350 220];
    
    S = fft(s);
    S = S/Fe;
    N = length(S);
    S = S([floor(N/2)+1:end, 1:floor(N/2)]);
    Df = Fe/length(t);
    f = (-N+floor(N/2):floor(N/2)-1) * Df;
    
    fig = figure;
    plot(f, abs(S).^2);
    xlabel('\nu');
    ylabel('|S(\nu)|^2');
    yticks(0);
    % ylim([-pi 0] + 0.*[-1 1]);
    xlim([0 17]);
    % ylim([0 1]);
    fig.Position(3:4) = [350 220];
end


S0 = nan(0, length(S));

for k = 1:200
    f = randn(size(t));
    
    s = F0*syst.response(f, dt, 1);
    
    S = fft(s);
    S = S/Fe;
    N = length(S);
    S = S([floor(N/2)+1:end, 1:floor(N/2)]);
    
    S0(end+1, :) = abs(S).^2;
end

S = mean(S0, 1);
Df = Fe/length(t);
f = (-N+floor(N/2):floor(N/2)-1) * Df;

fig = figure;
plot(f, S);
xlabel('\nu');
ylabel('|S(\nu)|^2');
yticks([0]);
% ylim([-pi 0] + 0.*[-1 1]);
xlim([0 17]);
% ylim([0 1]);
fig.Position(3:4) = [480 320];

%%

f0 = 10;
A = 1;
T = [0, 10000];
Fe = 500;
Fe2 = 154;

lambd = 2i*pi*f0;
t = T(1):1/Fe:T(2);
s = imag(A*exp(lambd*t));
t2 = T(1):1/Fe2:T(2);
s2 = imag(A*exp(lambd*t2));
T2 = [0.2, 0.8];
s2(~(t2>=T2(1) & t2<=T2(2))) = nan;

fig = figure;
plot(t, s);
hold on
scatter(t2, s2, 'r');
xline(T2(1), '--', 'LineWidth', 1);
xline(T2(2), '--', 'LineWidth', 1);
xlabel('t');
ylabel('s');
xlim([0 1]);
ylim(1.1*[-1 1]);
fig.Position(3:4) = [480 320];
xticks([]);
yticks([]);

