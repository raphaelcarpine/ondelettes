function BoundsQMenu(TLim, TLimRidge, ct, ctRidge, cf, MotherWavelet)
%BOUNDSQMENU Summary of this function goes here
%   Detailed explanation goes here


dlgtitle = 'Bounds for Q';
prompt = {'Frequency [Hz]:', 'Damping [%]:', 'Closest frequency(ies) [Hz]:'};
dims = [1 38];
answer = inputdlg(prompt, dlgtitle, dims);
try
    f = str2double(answer{1});
    zeta = str2double(answer{2})/100;
    f2 = str2num(answer{3});
catch
end

%%
if size(TLim, 1) > 1
    TLim = transpose(TLim);
end
if size(TLimRidge, 1) > 1
    TLimRidge = transpose(TLimRidge);
end
if size(ct, 1) > 1
    ct = transpose(ct);
end
if size(ctRidge, 1) > 1
    ctRidge = transpose(ctRidge);
end

%%

f2 = [f2, -f];  % cos(xt) = 1/2(exp(iwt) + exp(-iwt))


[Df, If2] = min(abs(f2 - f));
f2 = f2(If2);

[Qmin, Qmax, Qz] = getBoundsQ2(f, Df, zeta, TLim, TLimRidge, ct, ctRidge, cf, MotherWavelet);

%%

TitleMsg = 'Choice of Q';

message = {};
message{end+1} = 'Bounds for Q (Qmin < Q < Qmax, Q_z):         ';
message{end+1} = sprintf('Qmin: %.2f', Qmin);
message{end+1} = sprintf('Qmax: %.2f', Qmax);
message{end+1} = sprintf('Qz: %.2f', Qz);

message{end+1} = '';
message{end+1} = 'Infos:';
message{end+1} = sprintf('freq: %.2f Hz', f);
message{end+1} = sprintf('closest freq: %.2f Hz', f2);
message{end+1} = sprintf('damping: %.2f %%', 100*zeta);
message{end+1} = ['cf: ', num2str(cf)];
message{end+1} = ['time limits: ', num2str(TLim), ' s'];
message{end+1} = ['ct: ', num2str(ct)];
message{end+1} = ['time limits for ridge: ', num2str(TLimRidge), ' s'];
message{end+1} = ['ct for ridge: ', num2str(ctRidge)];
message{end+1} = ['mother wavelet: ', MotherWavelet];




msgbox(message, TitleMsg, 'help');


end

