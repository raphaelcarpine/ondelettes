load('pont sens/simulation elements finis/influence masse ajoutee/results');

%% valeurs pour mu <= 0.6

I = Array_mu <= 0.6;
Array_mu = Array_mu(I);
Array_Fn = Array_Fn(I);
Array_Fn_t = Array_Fn_t(I);
Array_Lw = Array_Lw(I);


%% tri

[Array_mu, I] = sort(Array_mu);
Array_Fn = Array_Fn(I);
Array_Fn_t = Array_Fn_t(I);
Array_Lw = Array_Lw(I);


%% affichage

error_freq = Array_Fn_t./Array_Fn .* sqrt(1 + Array_mu) - 1;
error_freq = abs(error_freq);

ind1 = Array_Lw == 1;
ind185 = Array_Lw == 18.5;

fig = figure;
plot(Array_mu(ind1), 100*error_freq(ind1), '-X');%, 'Linewidth', 1);
hold on
plot(Array_mu(ind185), 100*error_freq(ind185), '-+');%, 'Linewidth', 1);
xlim([0, 0.65]);
% plot([0, 1], [0, 0], '--', 'Color', 0.5*[1 1 1]);
xlabel('$\mu_t / \mu_b$', 'FontSize', 13,'Interpreter', 'latex');
ylabel('Error on frequency ratio [%]');
pos = get(fig, 'Position');
set(fig, 'Position', [pos(1:2), 400, 300]);
legend({'d_t = 1', 'd_t = 18.5'}, 'Location', 'NW');




