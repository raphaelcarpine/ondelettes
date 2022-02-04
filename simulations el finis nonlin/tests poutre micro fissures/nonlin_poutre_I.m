S = 1;
h = 1;
E1 = 1;
E2 = 0.5*E1;
sig0 = 1;

c0 = 2*sig0/(h*E1);


C = linspace(0, 3*c0 + 1*(c0==0), 10000);

x0 = ((1-E2/E1)*sig0 + C*h*E2)./(C*(E1+E2));

figure;
plot(C, x0);
hold on
% plot(C, x0bis);
xline(c0, '--');
yline(h/2, '--');
ylim([0, 1.5*h]);

%% figure courbure-moment

M0 = C*E1*h^2/2*S;
M1 = S *( C.*x0.^2*E1 + C.*(h-x0).^2*E2 + (h-x0)*(1-E2/E1)*sig0);
M1bis = S *( C*h^2*E1*E2/(E1+E2) + h*(E1-E2)/(E1+E2)*sig0);
M = M0.*(C <= c0) + M1.*(C > c0);
Mbis = M0.*(C <= c0) + M1bis.*(C > c0);

fig = figure;
plot([flip(-C), C], [flip(-M), M], 'LineWidth', 2);
hold on
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot([flip(-C), C], [flip(-M0), M0], '--', 'LineWidth', 2);
plot([flip(-C), C], [flip(-Mbis), Mbis], '-.', 'LineWidth', 2);
% plot( C, M2, '--');
% plot(flip(-C), flip(-M2), '--');
% xline(c0, '--');
plot(-c0*[1 1], [ax.YLim(1), -c0*E1*h^2/2*S], ':k');
plot(c0*[1 1], [ax.YLim(1), c0*E1*h^2/2*S], ':k');
xline(0, '-');
yline(0, '-');

xlim(max(C)*[-1 1]);

xticks([-c0, 0, c0]);
xticklabels({'-c_0', '0', 'c_0'});
yticks(0);
xlabel('Courbure')
ylabel('Moment')

fig.Position(3:4) = [380 300];

set(fig,'renderer','Painters');
% set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')


%% figure deformation-contrainte

sigma0 = 0.35;

fig = figure;
fig.Position(3:4) = [380 300];
plot([-1 sigma0 1], [-1 sigma0 sigma0+0.5*(1-sigma0)], 'LineWidth', 2);
hold on
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot([-1 1], [-1 1], '--', 'LineWidth', 2);
plot([-1 sigma0], [sigma0 sigma0], ':k');
fig.Position(3:4) = [380 300];
xlim([-1 1]);
ylim(1.1*[-1 1]);
xline(0, '-');
yline(0, '-');

xlabel('Déformation')
ylabel('Contrainte')
% xlabel('\epsilon', 'FontSize', 15)
% ylabel('\sigma', 'FontSize', 15)

xticks([0]);
yticks([0 sigma0]);
yticklabels({'0', '\sigma_0'});

text((-1+sigma0)/2, (-1+sigma0)/2, 'E', 'VerticalAlignment', 'top');
text((1+sigma0)/2, sigma0+0.5*(1-sigma0)/2, 'E''', 'VerticalAlignment', 'top');













