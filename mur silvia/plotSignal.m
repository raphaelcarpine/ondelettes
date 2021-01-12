P = 7;

%%
[t, X] = getData(P, 0);

figure;
plot(t, X);
xlabel('Time [s]');
ylabel('Acceleration [m/s²]');
xlim([5 65]);
l = [char(ones(9, 1) * 'CH'), num2str((1:9)')];
legend(l, 'NumColumns', 3, 'Location', 'best');