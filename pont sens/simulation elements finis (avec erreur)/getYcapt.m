function Ycapt = getYcapt(Ytot, pos_capt, dx)
%GETYCAPT Interpolation de la position des capteurs
%   Detailed explanation goes here

N = size(Ytot, 1);
Nt = size(Ytot, 2);
Ycapt = nan(length(pos_capt), Nt);
for it = 1:Nt
    Ycapt(:, it) = transpose(...
        Ytot( max(min( floor(pos_capt/dx) + 1, N), 1), it)' .* (1 - pos_capt/dx + floor(pos_capt/dx))...
        + Ytot( max(min( floor(pos_capt/dx) + 2, N), 1), it)' .* (pos_capt/dx - floor(pos_capt/dx)));
end
end

