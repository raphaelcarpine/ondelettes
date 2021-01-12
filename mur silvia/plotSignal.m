P = 6;

%%
[t, X] = getData(P, 0);

figure('Name', sprintf('Acceleration data, P%u', P));
plot(t, X);
xlabel('Time [s]');
ylabel('Acceleration [m/s²]');
% xlim([5 65]);
l = [char(ones(9, 1) * 'CH'), num2str((1:9)')];
legend(l, 'NumColumns', 3, 'Location', 'best');