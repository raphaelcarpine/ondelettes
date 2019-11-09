w2 = 1;
zeta = [0.01, 0.02, 0.04, 0.06];
legendZeta = cell(1, length(zeta));
for k = 1:length(zeta)
    z = zeta(k);
    legendZeta{k} = ['\zeta_2 = ', num2str(100*z), '%'];
end


w = linspace(0.9, 1.1, 1000);



fig = figure;
ax = axes(fig);
hold(ax, 'on');






for z = zeta
    H = @(p) p.^2 .* (2*z*w2*p + w2^2) ./ (p.^2 + 2*z*w2*p + w2^2);
    
    P = @(w) 1/2 * real(-1i*w .* H(1i*w));
    
    plot(ax, w, P(w));
    
    xlabel(ax, '\omega', 'FontSize', 15);
    ylabel(ax, 'P', 'FontSize', 15);
end

axis(ax, 'tight');

set(ax,'xtick',[0.9*w2, w2, 1.1*w2]);
set(ax,'xticklabels',{'0.9\times\omega_2', '\omega_2', '1.1\times\omega_2'});
set(ax,'ytick',[]);

legend(ax, legendZeta, 'Location', 'northwest', 'FontSize', 13);
legend(ax, 'boxoff')





