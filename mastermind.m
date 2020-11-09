n_c = 8; % nb de couleurs
n_p = 4; % nb de pions

couleurs = {'blanc', 'orange', 'bleu', 'rouge', 'vert', 'jaune', 'mauve', 'magenta'};

combinaisons_init = nan(n_c^n_p, n_p);
for k = 0:n_c^n_p-1
    combinaison_string = dec2base(k, n_c, n_p);
    for i = 1:n_p
        combinaisons_init(k+1, i) = str2num(combinaison_string(i));
    end
end



%% partie

% combinaisons = combinaisons_init;
% for k_tour = 1:12
%     fprintf('~~~~ tour %d ~~~~\n', k_tour);
%     [comb, esp] = meilleureComb(combinaisons_init, combinaisons, n_c, n_p);
%     fprintf('coup à jouer : (%d combinaisons restantes, %.1f%% d''elimination)\n',...
%         [size(combinaisons, 1), esp*100]);
%     disp(getCouleurs(comb, couleurs));
%
%     nbr = input('nombre de plots rouges : ');
%     nbb = input('nombre de plots blancs : ');
%     nb1 = nbr;
%     nb2 = nbr + nbb;
%
%     combinaisons = getNewCombinaisons(combinaisons, comb, nb2, nb1, n_c);
%
%     if size(combinaisons, 1) == 1
%         disp(' ');
%         disp('trouvé !')
%         disp(getCouleurs(combinaisons, couleurs));
%         return
%     end
%     if size(combinaisons, 1) == 0
%         disp(' ');
%         disp('hmmmmmm ? ? ?');
%         return
%     end
% end
%
% disp(' ');
% disp('pas trouvé :(')
% for k_comb = 1:size(combinaisons, 1)
%     disp(getCouleurs(combinaisons(k_comb, :), couleurs));
% end
% return


%% test

% comb_secrete = randi(n_c, 1, n_p) -1;
% disp('combinaison secrete :')
% disp(comb_secrete);
% disp(' ');
%
% combinaisons = combinaisons_init;
% for k_tour = 1:12
%     fprintf('~~~~ tour %d ~~~~\n', k_tour);
%     [comb, esp] = meilleureComb(combinaisons_init, combinaisons, n_c, n_p);
%     fprintf('coup à jouer : (%d combinaisons restantes, %.1f%% d''elimination)\n',...
%         [size(combinaisons, 1), esp*100]);
%     disp(comb);
%
%     nb1 = nbCouleurPlacement(comb_secrete, comb);
%     nb2 = nbCouleur(comb_secrete, comb, n_c);
%     fprintf('resultat : %d rouges, %d blancs\n\n', [nb1, nb2-nb1]);
%
%     combinaisons = getNewCombinaisons(combinaisons, comb, nb2, nb1, n_c);
%
%     if size(combinaisons, 1) == 1
%         disp('trouvé !')
%         disp(combinaisons);
%         return
%     end
% end


%% test 2

nb_tours = [];
proportion_tours = zeros(1, 8);
figure;
bar_graph = bar(proportion_tours);
ylim([0, 100]);
ytickformat('percentage');
txt_tour = text(1, 90, sprintf('nombre d''essais : %d', 0));
txts_prop = text(1:length(proportion_tours), proportion_tours, num2str(proportion_tours'),...
    'vert','bottom','horiz','center');

k_test = 1;
while true
    
    comb_secrete = randi(n_c, 1, n_p) -1;
    combinaisons = combinaisons_init;
    for k_tour = 1:100
        [comb, esp] = meilleureComb(combinaisons_init, combinaisons, n_c, n_p);
        
        nb1 = nbCouleurPlacement(comb_secrete, comb);
        nb2 = nbCouleur(comb_secrete, comb, n_c);
        
        if nb1 == n_p
            break
        end
        
        combinaisons = getNewCombinaisons(combinaisons, comb, nb2, nb1, n_c);
        
        if size(combinaisons, 1) == 0
            warning('pas de solution trouvée');
            break
        end
    end
    
    %     nb_tours = [nb_tours, k_tour];
    %     fprintf('nombre de tours moyen : %.2f, max : %d, min %d (%d essais)\n',...
    %         [mean(nb_tours), max(nb_tours), min(nb_tours), length(nb_tours)]);
    
    proportion_tours(k_tour) = proportion_tours(k_tour) + 1;
    set(bar_graph, 'YData', proportion_tours/k_test*100);
    set(txt_tour, 'String', sprintf('nombre d''essais : %d', k_test));
    for k_prop = 1:length(proportion_tours)
        set(txts_prop(k_prop), 'String', sprintf('%.1f%%', proportion_tours(k_prop)/k_test*100), ...
            'Position', [k_prop, proportion_tours(k_prop)/k_test*100]);
    end
    drawnow;
    
    k_test = k_test + 1;
end







%%

function str = getCouleurs(comb, couleurs)
str = '';
for c = comb
    str = [str, couleurs{c+1}, ', '];
end
str = str(1:end-2);
end

function nb = nbCouleurPlacement(combinaisons, combinaison_test)
nb = (combinaisons - combinaison_test) == 0;
nb = sum(nb, 2);
end


function nb = nbCouleur(combinaisons, combinaison_test, n_c)
nb = zeros(size(combinaisons, 1), 1);

for couleur = 0:n_c-1
    nb_couleur_test = sum(combinaison_test == couleur, 2);
    if nb_couleur_test == 0
        continue
    end
    nb = nb + min(sum(combinaisons == couleur, 2), nb_couleur_test);
end
end


function combinaisons = getNewCombinaisons(combinaisons, comb_test, nb_couleurs, nb_placements, n_c)
nb1 = nbCouleurPlacement(combinaisons, comb_test);
nb2 = nbCouleur(combinaisons, comb_test, n_c);
combinaisons = combinaisons((nb1==nb_placements)&(nb2==nb_couleurs), :);
end

function [comb, esperance] = meilleureComb(combinaisons_init, combinaisons, n_c, n_p)
esperance_elimination = nan(n_c^n_p, 1);
for k_test = 1:n_c^n_p
    comb_test = combinaisons_init(k_test, :);
    
    nb1 = nbCouleurPlacement(combinaisons, comb_test);
    nb2 = nbCouleur(combinaisons, comb_test, n_c);
    nb = n_p*nb1 + nb2;
    nb = sort(nb);
    
    grps_identiques = [1];
    for k_nb = 2:length(nb)
        if nb(k_nb-1) == nb(k_nb)
            grps_identiques(end) = grps_identiques(end) + 1;
        else
            grps_identiques = [grps_identiques, 1];
        end
    end
    
    esperance_elimination(k_test) = ...
        sum(grps_identiques .* (length(nb)-grps_identiques));
end

esperance = max(esperance_elimination);
K_comb = esperance_elimination == esperance;
K_comb_ind = [];
for k_comb = 1:length(K_comb)
    if K_comb(k_comb)
        K_comb_ind = [K_comb_ind, k_comb];
    end
end
combs = combinaisons_init(K_comb, :);
ind_combs_in_combinaisons = false(size(combs, 1), 1);
for k_comb = 1:size(combs, 1)
    ind_combs_in_combinaisons(k_comb) = any( all( combinaisons == combs(k_comb, :), 2));
end
if any(ind_combs_in_combinaisons)
    combs = combs(ind_combs_in_combinaisons, :);
end

comb = combs(randi(size(combs, 1)), :);

esperance = esperance / length(nb)^2;
end











function [comb, esperance] = meilleureComb2(combinaisons_init, combinaisons, n_c, n_p)
esperance_elimination = nan(n_c^n_p, 1);
for k_test = 1:n_c^n_p
    comb_test = combinaisons_init(k_test, :);
    
    nb1 = nbCouleurPlacement(combinaisons, comb_test);
    nb2 = nbCouleur(combinaisons, comb_test, n_c);
    nb = n_p*nb1 + nb2;
    [nb, Icomb] = sort(nb);
    
    nb_identiques = ones(1, length(nb));
    for k_nb = 2:length(nb)
        if nb(k_nb-1) == nb(k_nb)
            nb_identiques(k_nb) = nb_identiques(k_nb-1) + 1;
        end
    end
    
    grps_identiques = [];
    k_nb = length(nb_identiques);
    while k_nb > 0
        nb_id = nb_identiques(k_nb);
        
        grps_identiques = [grps_identiques, nb_id];
        
        nb_identiques(k_nb-nb_id+1:k_nb) = nb_id * ones(nb_id, 1);
        k_nb = k_nb - nb_id;
    end
    
    esperance_elimination(k_test) = ...
        sum(grps_identiques .* (length(nb)-grps_identiques))/length(nb)^2;
end

[esperance, k_comb] = max(esperance_elimination);
comb = combinaisons_init(k_comb, :);
end