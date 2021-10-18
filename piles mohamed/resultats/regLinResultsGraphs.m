dataFolder = 'piles mohamed\resultats\reg lin';
saveFolder = 'piles mohamed\resultats\figures reg lin';

saveFigs = true;

pile = 10;
direction = 'x';
g = 8;
mode = 2;

%%

fileName = sprintf('pile%u_%c_%ug_mode%u_reglin.xlsx', [pile, int8(direction), g, mode]);
A = readmatrix(fullfile(dataFolder, fileName));

A(:, 4:5) = 100*A(:, 4:5); % pourcentage

% freq
prof = A(~isnan(A(:, 2)), 1);
f0 = A(~isnan(A(:, 2)), 2);
fslope = A(~isnan(A(:, 2)), 3);
fig = figure;
yyaxis left
plot(prof, f0, '--');
hold on
plot(prof, f0, '+', 'LineWIdth', 2, 'MarkerSize', 10);
ylabel('Frequency intercept [Hz]');
yyaxis right
plot(prof, fslope, '--');
hold on
plot(prof, fslope, '+', 'LineWIdth', 2, 'MarkerSize', 10);
ylabel('Frequency slope [Hz.s²/m]');
xlabel('Depth [cm]');
xlim(A([1 end], 1).' + 1*[-1 1]);
xticks(A(:, 1));

if saveFigs
    fileName = sprintf('pile%u_%c_%ug_mode%u_reglin_freq', [pile, int8(direction), g, mode]);
    saveas(fig, fullfile(saveFolder, fileName), 'fig');
    set(fig,'renderer','Painters');
    saveas(fig, fullfile(saveFolder, fileName), 'epsc');
    saveas(fig, fullfile(saveFolder, fileName), 'png');
    fprintf('''%s'' saved\n', fileName);
end

% damp
prof = A(~isnan(A(:, 2)), 1);
damp0 = A(~isnan(A(:, 2)), 4);
dampslope = A(~isnan(A(:, 2)), 5);
fig = figure;
yyaxis left
plot(prof, damp0, '--');
hold on
plot(prof, damp0, '+', 'LineWIdth', 2, 'MarkerSize', 10);
ylabel('Damping intercept [%]');
yyaxis right
plot(prof, dampslope, '--');
hold on
plot(prof, dampslope, '+', 'LineWIdth', 2, 'MarkerSize', 10);
ylabel('Damping slope [%.s²/m]');
xlabel('Depth [cm]');
xlim(A([1 end], 1).' + 1*[-1 1]);
xticks(A(:, 1));

if saveFigs
    fileName = sprintf('pile%u_%c_%ug_mode%u_reglin_damp', [pile, int8(direction), g, mode]);
    saveas(fig, fullfile(saveFolder, fileName), 'fig');
    set(fig,'renderer','Painters');
    saveas(fig, fullfile(saveFolder, fileName), 'epsc');
    saveas(fig, fullfile(saveFolder, fileName), 'png');
    fprintf('''%s'' saved\n', fileName);
end