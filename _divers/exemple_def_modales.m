pile = 1;
% acceleromètres (acc1 tout en haut, 2 juste en dessous etc.)
chNames = {'acc1x', 'acc1y', 'acc1z', 'acc2x', 'acc2y', 'acc3x', 'acc3y', 'acc4x', 'acc4y'}; % 'acc5x', 'acc5y'}
nameFig = '';
zSable = 10;
direction = 'x'; % direction choc, pour deformee 1D

deformee = 0.2*(randn(1, length(chNames)) + 1i*randn(1, length(chNames))); % deformee au hasard
deformee = deformee / sqrt(deformee * deformee.'); % normalisation dans C


[fctDefModale, fctDefModaleAnimation] = defModales(pile, chNames, nameFig, zSable);
fctDefModale1D = defModales1D(pile, chNames, direction, nameFig, zSable);

fctDefModale(real(deformee), 'deformee 3D'); % rouge: x, violet: y, bleu: z
fctDefModaleAnimation(deformee, 'deformee 3D animation');
fctDefModale1D(real(deformee), 'deformee 1D');