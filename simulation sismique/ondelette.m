


% data = xlsread('regular.xlsx');
% data = xlsread('irregular.xlsx');
data = xlsread('irregular2.xlsx');


%%


t = data(:,1);
t = t(~isnan(t));
t = t(1:end-1);

A0 = data(:,2);
A0 = A0(~isnan(A0));
A0 = A0(1:end-1);

tx = data(:,3);
x = data(:,4);
v = data(:,6);
a = data(:,8);

X = nan(size(t));
V = nan(size(t));
A = nan(size(t));

itx = 1;
for it=1:length(t)
    while itx < length(tx) && tx(itx+1) == tx(itx)
        itx = itx+1;
    end
    
    if tx(itx) ~= t(it)
        warning('');
    end
    
    X(it) = x(itx);
    V(it) = v(itx);
    A(it) = a(itx);
    
    itx = itx+1;
end

%% choix reponse libre

n0 = length(t); % debut de la reponse libre
while A0(n0-1) == 0
    n0 = n0-1;
end

% t = t(n0:end);
% A0 = A0(n0:end);
% X = X(n0:end);
% V = V(n0:end);
% A = A(n0:end);

t = t(1:n0);
A0 = A0(1:n0);
X = X(1:n0);
V = V(1:n0);
A = A(1:n0);

%% filtrage

X = X - mean(X);
V = V - mean(V);
A = A - mean(A);


% fX = fft(X);
% fV = fft(V);
% fA = fft(A);
% 
% f0 = 10; % fréquance de filtrage passe haut
% 
% dt0 = t(2)-t(1);
% fshanon = 1/(2*dt0);
% n0 = round(length(fX)/2 * f0/fshanon);
% 
% fX(1:n0) = zeros(n0, 1);
% fX(end-n0+1:end) = zeros(n0, 1);
% fX(floor(length(fX)/2):end) = zeros(floor(length(fX)/2)+2, 1);
% 
% X = ifft(fX);




%% plot

% X = ones(size(X));

% figure;
% pa0 = plot(t, A0);
% xlabel('t');
% ylabel('a0');
% 
% figure;
% px = plot(t, X);
% xlabel('t');
% ylabel('x');
% 
% figure;
% pv = plot(t, V);
% xlabel('t');
% ylabel('v');
% 
% figure;
% pa = plot(t, A);
% xlabel('t');
% ylabel('a');


%% ondelette


Q = 10;
MaxRidges = 100;
MaxParallelRidges = inf;
fmin = 2;
fmax = 6;
NbFreq = 100;
WvltFreq = linspace(fmin, fmax, NbFreq);
ctEdgeEffects = 3;


wvltA= WvltComp(t, A , WvltFreq, Q);
wvltA0= WvltComp(t, A0 , WvltFreq, Q);

wvltH = wvltA./wvltA0;


f = figure;
axH = axes(f);
pcolor(t, WvltFreq, log10(abs(wvltH)));
xlabel('time [T]')
ylabel('frequency [T]^{-1}')
shading flat
colormap(jet)
hold on
delta_t = ctEdgeEffects*Q./(2*pi*WvltFreq); % = a * Delta t_psi, en considerant que 2*pi*Delta f_psi * Delta t_psi = 1/2
plot(t(1)+delta_t,WvltFreq,'black','LineWidth',1.5,'LineStyle',':') % limite effets de bord gauche
plot(t(end)-delta_t,WvltFreq,'black','LineWidth',1.5,'LineStyle',':') % limite effets de bord droite
hold off
title('H');


f = figure;
axA0 = axes(f);
pcolor(t, WvltFreq, log10(abs(wvltA0)));
xlabel('time [T]')
ylabel('frequency [T]^{-1}')
shading flat
colormap(jet)
hold on
delta_t = ctEdgeEffects*Q./(2*pi*WvltFreq); % = a * Delta t_psi, en considerant que 2*pi*Delta f_psi * Delta t_psi = 1/2
plot(t(1)+delta_t,WvltFreq,'black','LineWidth',1.5,'LineStyle',':') % limite effets de bord gauche
plot(t(end)-delta_t,WvltFreq,'black','LineWidth',1.5,'LineStyle',':') % limite effets de bord droite
hold off
title('a0');


f = figure;
axA = axes(f);
pcolor(t, WvltFreq, log10(abs(wvltA)));
xlabel('time [T]')
ylabel('frequency [T]^{-1}')
shading flat
colormap(jet)
hold on
delta_t = ctEdgeEffects*Q./(2*pi*WvltFreq); % = a * Delta t_psi, en considerant que 2*pi*Delta f_psi * Delta t_psi = 1/2
plot(t(1)+delta_t,WvltFreq,'black','LineWidth',1.5,'LineStyle',':') % limite effets de bord gauche
plot(t(end)-delta_t,WvltFreq,'black','LineWidth',1.5,'LineStyle',':') % limite effets de bord droite
hold off
title('a');





%% ridges

ridgesA0 = RidgeExtract(t, nan, Q, fmin, fmax, NbFreq, 'Wavelet', wvltA0, 'NbMaxRidges', MaxRidges,...
    'NbMaxParallelRidges', MaxParallelRidges, 'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects);

ridgesA = RidgeExtract(t, nan, Q, fmin, fmax, NbFreq, 'Wavelet', wvltA, 'NbMaxRidges', MaxRidges,...
    'NbMaxParallelRidges', MaxParallelRidges, 'ctLeft', ctEdgeEffects, 'ctRight', ctEdgeEffects);


f = figure;
ax = axes(f);
f = figure;
axAmpl = axes(f);
hold(ax, 'on');
hold(axAmpl, 'on');
hold(axA0, 'on');
hold(axA, 'on');
for k = 1:length(ridgesA0.time)
    plot(ridgesA0.time{k}, ridgesA0.freq{k}, 'b', 'Parent', ax);
    plot(ridgesA0.time{k}, abs(ridgesA0.val{k}), 'b', 'Parent', axAmpl);
    plot(ridgesA0.time{k}, ridgesA0.freq{k}, 'black', 'Parent', axA0);
end
for k = 1:length(ridgesA.time)
    plot(ridgesA.time{k}, ridgesA.freq{k}, 'r', 'Parent', ax);
    plot(ridgesA.time{k}, abs(ridgesA.val{k}), 'r', 'Parent', axAmpl);
    plot(ridgesA.time{k}, ridgesA.freq{k}, 'black', 'Parent', axA);
end
hold(ax, 'off');
hold(axAmpl, 'off');
hold(axA0, 'off');
hold(axA, 'off');
xlim(ax, [t(1), t(end)]);
ylim(ax, [fmin, fmax]);
xlim(axAmpl, [t(1), t(end)]);









% QtyX = {
%     'time'
%     'time'
%     'time'
%     };
% 
% QtyY = {
%     'freq'
%     'val'
%     'bandwidth'
%     };
% 
% FuncX = {
%     ''
%     ''
%     ''
%     };
% 
% FuncY = {
%     ''
%     'abs'
%     ''
%     };
% 
% 
% 
% 
% for k = 1:length(QtyX)
%     RidgeQtyPlot2(ridges, QtyX{k}, QtyY{k} , 'EvaluationFunctionX', FuncX{k}, 'EvaluationFunctionY', FuncY{k});
% end
% 
% 
% WaveletMenu('WaveletPlot', pa, 'fmin', fmin, 'fmax', fmax,...
%     'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);






