function [Tdj1, Tdj2, Tjour] = getData(jour)

if nargin == 0
    jour = 1;
end

fileFolder = 'C:\Users\carpine\Documents\projets\donnees poids lourds franziska\donnees jours semaine';
fileName = sprintf('%d.mat', jour);

load(fullfile(fileFolder, fileName));

%%

dj1 = datetime(Tjour.date(1)) - timeofday(Tjour.date(1)); % date
dj2 = dj1 + 0.5; % date midi

Tdj1 = Tjour(Tjour.date < dj2, :); % première demie-journée
Tdj2 = Tjour(Tjour.date >= dj2, :); % deuxième demie-journée

Tdj1.time = seconds(Tdj1.date - dj1);
Tdj2.time = seconds(Tdj2.date - dj2);

end