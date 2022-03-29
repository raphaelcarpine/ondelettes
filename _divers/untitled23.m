%% loi arrivee vehicules

t = linspace(-5, 35, 1000);
p = (t>=0) .* exp(-t/6);

fig = figure;
plot(t, p);
ylim([0, 1.15]);
xlabel('Temps [s]');
ylabel('Densité de probabilité');
yticks([]);
fig.Position(3:4) = [300 220];

%% loi masse

% conversion loi log-normale
sigma2_m_log = log((250/1250)^2 + 1);
mu_m_log = log(1250) - sigma2_m_log/2;

m = linspace(0, 2500, 1000);
p = exp(-(log(m)-mu_m_log).^2/(2*sigma2_m_log)) ./ m;

fig = figure;
plot(m, p);
ylim([0, 1.15*max(p)]);
xlabel('Masse [kg]');
ylabel('Densité de probabilité');
yticks([]);
fig.Position(3:4) = [300 220];

%% loi vitesse

v = linspace(0, 100, 1000);
p = exp(-(v-55).^2/(2*10^2));

fig = figure;
plot(v, p);
ylim([0, 1.15]);
xlabel('Vitesse [km/h]');
ylabel('Densité de probabilité');
yticks([]);
fig.Position(3:4) = [300 220];