clear all
close all

%%
Ftrain = 4.2;
AmplTrain = [1.2, 1, 0.2];
PhaseTrain = [0, pi, pi];

AmplH = [0.2, 1, 2];
PhaseH = [0.2, 1, 2];

Fpropre = [6.5, 14.75];
Zeta = [0.05, 0.05];
Ampl0 = [0.5, 0.2];
Phase0 = 2*pi*rand(1, 3);
modes = [1 3];

T = 10;
t1 = 2;
t2 = t1 + 10/Ftrain;
t = 0:0.001:T;

folderPath = 'pont sens\figures bruno\figs\';


%% excitation

f = [];
for k = 1:length(AmplTrain)
    f(k, :) = AmplTrain(k) .* cos((k-1)*Ftrain*2*pi*(t-t1) + PhaseTrain(k)) .* (t>=t1 & t<=t2);
    fig = plot2(t, f(k, :), 'F', [-1.2, 2.2], T, [t1, t2]);
    saveas(fig, [folderPath, 'train_', num2str(k-1), 'F.png']);
end

fig = plot2(t, sum(f, 1), 'F', [-1.2, 2.2], T, [t1, t2]);
saveas(fig, [folderPath, 'train_tot.png']);

%% réponse harmonique

xh = [];
for k = 1:length(AmplTrain)
    xh(k, :) = AmplTrain(k)*AmplH(k) * cos((k-1)*Ftrain*2*pi*(t-t1) + PhaseTrain(k)+PhaseH(k)) .* (t>=t1 & t<=t2);
    fig = plot2(t, xh(k, :), 'x', [-1.2, 1.2], T, [t1, t2]);
    saveas(fig, [folderPath, 'reponse_harmonique_', num2str(k-1), 'F.png']);
end

fig = plot2(t, sum(xh, 1), 'x', [-1.2, 2.2], T, [t1, t2]);
saveas(fig, [folderPath, 'reponse_harmonique_tot.png']);


%% réponse libre

xl1 = [];
xl2 = [];
for k = 1:length(Fpropre)
    xl1(k, :) = Ampl0(k) * cos(2*pi*Fpropre(k)*(t-t1) + Phase0(k)) .* exp(-Zeta(k)*2*pi*(t-t1)) .* (t>=t1);
    fig = plot2(t, xl1(k, :), 'x', [-1.2, 1.2], T, [t1, t2]);
    saveas(fig, [folderPath, 'reponse_libre_mode', num2str(k), '_t1.png']);
    xl2(k, :) = - Ampl0(k) * cos(2*pi*Fpropre(k)*(t-t2) + Phase0(k)) .* exp(-Zeta(k)*2*pi*(t-t2)) .* (t>=t2);
    fig = plot2(t, xl2(k, :), 'x', [-1.2, 1.2], T, [t1, t2]);
    saveas(fig, [folderPath, 'reponse_libre_mode', num2str(k), '_t2.png']);
    fig = plot2(t, xl1(k, :)+xl2(k, :), 'x', [-1.2, 1.2], T, [t1, t2]);
    saveas(fig, [folderPath, 'reponse_libre_mode', num2str(k), '_t1t2.png']);
end






















function fig = plot2(x, y, yLabel, yLim, T, t1t2)
fig = figure;
plot(x, y, 'LineWidth', 1);
xlim([0, T]);
ylim(yLim);
xticks(t1t2);
xticklabels({'t_1', 't_2'});
yticks([])
set(gca, 'box', 'off');
ylabel(yLabel);
end