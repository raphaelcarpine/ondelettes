%% system
m = 1;
k1 = (10*2*pi)^2;
k2 = k1/2;
x0 = 1;
zeta = 0.005;
zeta = 0.01;
c = 2*zeta*sqrt(k1*m);


g = @(x) (k2-k1)*(x-x0).*(x>x0);
% g = @(x) (k2-k1)*(x-x0).*(x>x0) + (k2-k1)*(x+x0).*(x<-x0);
% g = @(x) 0;

if 0
    X = linspace(-3*x0, 3*x0, 1000);
    Kx = nan(size(X));
    Gx = nan(size(X));
    for ix = 1:length(X)
        Kx(ix) = k1*X(ix) + g(X(ix));
        Gx(ix) = g(X(ix));
    end
    figure;
%     xline(0);
%     hold on
%     yline(0);
    plot(X, Kx, 'LineWidth', 2);
%     hold on
%     plot(X, Gx);
    hold on
    plot([-x0, -x0], [get(gca, 'YLim') * [1;0], -k1*x0], '--', 'Color', [0 0 0]);
    plot([x0, x0], [get(gca, 'YLim') * [1;0], k1*x0], '--', 'Color', [0 0 0]);
    xlabel('Displacement');
    ylabel('Force');
    text(0, 0, {'k  ', ''}, 'FontSize', 14, 'HorizontalAlignment', 'right');
    xkp = 2*x0;
    text(xkp, k1*xkp+g(xkp), {'k'' ', ''}, 'FontSize', 14, 'HorizontalAlignment', 'right');
    text(-xkp, -k1*xkp+g(-xkp), {'k'' ', ''}, 'FontSize', 14, 'HorizontalAlignment', 'right');
    xticks([-x0, 0, x0]);
    xticklabels({'-x_0', '0', 'x_0'});
    yticks([]);
    xlim(3*x0*[-1 1]);
%     grid on
end


% time
T = 10000;
% T = 30;
T0 = 5/(c/(2*m));
fe = 1000;

dt = 1/fe;
T0 = dt * ceil(T0/dt);
t = -T0:dt:T;
kt0 = sum(t < 0) + 1;

%% excitation

% dirac
f = zeros(size(t));
f(kt0) = 1000 / dt;

% noise
f0 = 50;
f = f0/sqrt(dt) * randn(size(t));

% noise + static
f0 = 5;
f = f0/sqrt(dt) * randn(size(t)) + 2*k1*sin(2*pi*t/100);

% Aref
Aref = f0*sqrt(dt)/(2*m*sqrt(k1/m)*sqrt(zeta*sqrt(k1/m)*dt))

% sigma F
sigmaF = sqrt(f0^2*dt^2*pi/(m*k1/m*(1-zeta^2)*dt))


%% ode
xi = 0;
vi = 0;

Xi = [xi; vi];

D = @(t, X) [X(2); 1/m*(-k1*X(1) -g(X(1)) - c*X(2) - f(floor((T0+t)/dt)+1))];

X = RK4(D, t, Xi);
x = X(1, :);

% delete begining
t = t(kt0:end);
x = x(kt0:end);

% resampling
fe2 = 100;
Nfe = floor(fe/fe2);
t = t(1:Nfe:end);
x = x(1:Nfe:end);

%%

% plot
figure;
plt = plot(t, x);
xlabel('Time [s]');
ylabel('Displacement [m]');

% % plot hilbert
% figure('Name', 'Hilbert');
% hilb = hilbert(x);
% Dhilb = angle(hilb(2:end)./hilb(1:end-1));
% Dhilb = [Dhilb(1), (Dhilb(1:end-1)+Dhilb(2:end))/2, Dhilb(end)];
% plot(abs(hilb), Dhilb/(2*pi*dt));
% xlabel('Amplitude [m]');
% ylabel('Frequency [Hz]');
% 
% % plot hilbert 2
% figure('Name', 'Hilbert 2');
% P = round(2*pi/(sqrt(k1/m)*dt));
% hilb = hilbert(x);
% Dhilb = angle(hilb(1+P:end)./hilb(1:end-P));
% plot(abs(hilb(1+round(P/2):round(P/2)+length(Dhilb))), Dhilb/(2*pi*P*dt));
% xlabel('Amplitude [m]');
% ylabel('Frequency [Hz]');
% 
% 
% % plot CWT slope
% fmin = 5;
% fmax = 15;
% Q = 2;
% MotherWavelet = 'morlet';
% freqs = linspace(fmin, fmax, 100);
% CWT = WvltComp(t, x, freqs, Q, 'MotherWavelet', MotherWavelet);
% for slopeTime = [0.1 0.5 1 2]
%     ridge = SingleRidgeExtract(t, freqs, CWT, MotherWavelet, Q, 3, 'slope', slopeTime);
%     figure('Name', sprintf('CWT slope penalization %f', slopeTime));
%     plot(abs(ridge.val), ridge.freq);
%     xlabel('Amplitude [m]');
%     ylabel('Frequency [Hz]');
% end


%%
% wavelet menu

fmin = 3;
fmax = 15;
Q = sqrt(pi/zeta)/3;
MotherWavelet = 'morlet';
MaxParallelRidges = 1;
MaxSlopeRidge = inf;
SignalUnit = 'm';
SquaredSignalUnit = 'm²';


WaveletMenu('WaveletPlot', plt, 'fmin', fmin, 'fmax', fmax, 'Q', Q, 'MotherWavelet', MotherWavelet,...
    'MaxParallelRidges', MaxParallelRidges, 'MaxSlopeRidge', MaxSlopeRidge,...
    'SignalUnit', SignalUnit, 'SquaredSignalUnit', SquaredSignalUnit);


%% freq

f1 = @(A) 1/(2*pi) * (sqrt(k1/m) - (A > x0) *(k1-k2)/(pi*sqrt(m*k1)) .* ( ...
    acos(x0./A) - (x0./A) .* sqrt(1-(x0./A).^2)));

f2 = @(A) 1/(2*pi) * ( (A <= x0) *sqrt(k1/m) + (A > x0) * (pi/2) ./ ...
    (sqrt(m/k1)*asin(x0./sqrt(A.^2+(k2-k1)/k1*(A-x0).^2)) + sqrt(m/k2)*acos(k1*x0./(k2*A+(k1-k2)*x0))));


if false
    A = linspace(0, 10*x0, 1000);
    figure;
    plot(A, f1(A));
    hold on
    plot(A, f2(A));
end


Am = linspace(0, 2.5, 300);
Aasym = nan(size(Am));
fasym = nan(size(Am));
for k = 1:length(Am)
    [Aasym(k), fasym(k)] = f2asym(Am(k), m, k1, k2, x0);
end


function [A, f] = f2asym(Am, m, k1, k2, x0)
if Am <= x0
    A = Am;
    f = sqrt(k1/m)/(2*pi);
    return
end

d = Am - x0;
c1 = k2;
c2 = k1*2*x0;
c3 = -k1*(2*x0*d+d^2);
D = (-c2+sqrt(c2^2-4*c1*c3))/(2*c1);
T = sqrt(m/k1)*pi/2;
T = T + sqrt(m/k1)*atan(x0/sqrt(Am^2-x0^2));
T = T + sqrt(m/k2)*(atan((k1/k2*x0+D)/sqrt(k1/k2*2*x0*D+(k1/k2*x0)^2+D^2-(k1/k2*x0+D)^2)) -...
    atan((k1/k2*x0)/sqrt(k1/k2*2*x0*D+(k1/k2*x0)^2+D^2-(k1/k2*x0)^2)));
T = 2*T;
f = 1/T;
A = x0 + (d+D)/2;

end
















