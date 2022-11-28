function projMode = projectionMode(defMode, defsParasites)
%PROJECTION Permet de trouver la projection orthogonale idéale permettant
%d'isoler un mode, en éliminant les modes parasites.
%   n : nombre de capteurs, i.e. dimension des déformées
%   m : nombre de modes à éliminer
%   defMode : vecteur n*1 (ou 1*n), déformée modale du mode à isoler
%   defsParasites : matrice n*m (ou m*n), contenant les déformées des modes
%   à éliminer côte à côte (ou les unes sur les autres)
%   projMode : vecteur n*1 (ou 1*n), vecteur permettant de faire la
%   projection orthogonale avec projMode.'*signal (ou projMode*signal.')

%% mise en forme des vecteurs

% vérification orientation defMode
if iscolumn(defMode)
    retournement = false; % cas vecteurs colone
elseif isrow(defMode)
    retournement = true; % cas vecteurs ligne
    defMode = defMode.';
    defsParasites = defsParasites.';
else
    error('defMode doit être un vecteur');
end

% tailles
n = size(defMode, 1); % dimension des déformées
m = size(defsParasites, 2); % nombre de modes à éliminer

if size(defsParasites, 1) ~= n % verification taille déformées parasites
    error('la dimension de defsParasites doit être n*m (ou m*n)');
end

if m > n-1
    error('le nombre de déformées parasites m ne peut pas excéder n-1');
end

% normalisation dans le plan complexe
defMode = defMode / sqrt(defMode.' * defMode);
for k = 1:m
    defsParasites(:, k) = defsParasites(:, k) / sqrt(defsParasites(:, k).' * defsParasites(:, k));
end

% partie réelle
defMode = real(defMode);
defsParasites = real(defsParasites);


%% construction des projections possibles

while true
    M = [defsParasites'; randn(n-m, n)]; % matrice des projections
    
    Psi = M\[zeros(m, n-m); eye(n-m)]; % ensemble de projetées orthogonales
    A0 = (Psi.'*Psi)\(Psi.'*defMode); % vecteur maximisant sa projection avec defMode (de norme 1), dans la base P
    projMode = Psi*A0; % vecteur maximisant sa projection avec defMode
    projMode = projMode / sqrt(projMode.'*projMode); % projetée maximisant proj0.'*defMode pour proj0.'*proj0 = 1
    
    if ~any(isnan(projMode)) % pb de conditionnement, possible à cause du tirage aléatoire
        break
    end
end

%% sortie

if retournement % remise en forme, dans le cas d'une entrée en vecteurs lignes
    projMode = projMode.';
end


end

