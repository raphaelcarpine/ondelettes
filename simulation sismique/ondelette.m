


% data = xlsread('regular.xlsx');
data = xlsread('irregular.xlsx');


%%


t = data(:,1);
t = t(~isnan(t));
t = t(1:end-1);

x0 = data(:,2);
x0 = x0(~isnan(x0));
x0 = x0(1:end-1);

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
while x0(n0-1) == 0
    n0 = n0-1;
end

% t = t(n0:end);
% x0 = x0(n0:end);
% X = X(n0:end);
% V = V(n0:end);
% A = A(n0:end);

t = t(1:n0);
x0 = x0(1:n0);
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

figure;
px0 = plot(t, x0);
xlabel('t');
ylabel('x0');

figure;
px = plot(t, X);
xlabel('t');
ylabel('x');

% figure;
% pv = plot(t, V);
% xlabel('t');
% ylabel('v');

% figure;
% pa = plot(t, A);
% xlabel('t');
% ylabel('a');


%% ondelette


Q = 5;
MaxParallelRidges = 1;
fmin = 4;
fmax = 6;
NbFreq = 300;


WaveletMenu('WaveletPlot', px, 'fmin', fmin, 'fmax', fmax,...
    'NbFreq', NbFreq, 'Q', Q, 'MaxParallelRidges', MaxParallelRidges);






