load('mur silvia\modesOndelette\allModalQuantities.mat');
P = [0, 6, 7];

%% freqs
disp('~~~~~~~~~~~~ Frequencies ~~~~~~~~~~~~');
for p = 1:3
    fprintf('~~~~~~ P%u ~~~~~~\n', P(p));
    for k = 1:length(AllModalQuantities.freqs{p})
        f = AllModalQuantities.freqs{p}{k};
        fprintf('\nk = %u : %.2fHz', [k, mean(f)]);
        if length(f) > 1
            fprintf(' +- %.2fHz (%.2f%%)', [std(f), 100*std(f)/mean(f)]);
        end
    end
    disp(' ');
end

%% damps
disp(' ');
disp(' ');
disp('~~~~~~~~~~~~ Damping ratios ~~~~~~~~~~~~');
for p = 1:3
    fprintf('~~~~~~ P%u ~~~~~~', P(p));
    for k = 1:length(AllModalQuantities.damps{p})
        d = AllModalQuantities.damps{p}{k};
        fprintf('\nk = %u : %.2f%%', [k, 100*mean(d)]);
        if length(d) > 1
            fprintf(' +- %.2f%% (%.1f%%)', [100*std(d), 100*std(d)/mean(d)]);
        end
    end
    disp(' ');
    disp(' ');
end
