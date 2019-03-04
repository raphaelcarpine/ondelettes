function WvltOut=raindrop(WvltIn)

modu = abs(WvltIn); % On ne travaille qu'en module
%%
[ny,nx]=size(modu); % ny = nb. freq, nx = nb. instants

%%
L_pos = 1:(nx*ny); % liste des indices des elements de modu
%%

L_pos_do = L_pos -1; % liste des indices des elements tout juste en dessous de L_pos (freq. tout juste inferieure)
L_pos_up = L_pos + 1; % liste des indices des elements tout juste au dessus de L_pos (freq. tout juste superieure)

L_pos_do(1:ny:(nx*ny))=1:ny:(nx*ny); % correction en bord inferieur de modu : il ne faut pas sortir de la colonne consideree
L_pos_up(ny:ny:(nx*ny))=ny:ny:(nx*ny);% correction en bord superieur de modu : il ne faut pas sortir de la colonne consideree


mod0 = modu(L_pos); % liste des elements de modu

modup = modu(L_pos_up); % liste des elements de modu de la frequence tout juste sup.
moddo = modu(L_pos_do); % liste des elements de modu de la frequence tout juste inf.

nomvt = (((modup-mod0)<0) & ((moddo-mod0)<0)); % Si modu plus grand que son element aux freq. tout juste sup. et inf. alors c'est une max. local


% L final
L_f = L_pos(nomvt); % On recupere la liste des indices des elements identifies comme max. locaux

%% mesu

WvltOut = zeros([ny,nx]); % Matrice meme taille que wvlt_in, zero partout

WvltOut(L_f) = WvltIn(L_f); % zero partout SAUF aux max. locaux (valeur de wvlt_in dans ce cas)
