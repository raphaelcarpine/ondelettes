w0 = 1; % \omega_0
z = 0.05; % \zeta
U = 1;
phi = 0;



l = z*w0; % \lambda
wa = w0*sqrt(1-z^2); % \omega_a


Fu = @(w) 1/2 * (exp(1i*phi)./(l + 1i*(w-wa)) + exp(-1i*phi)./(l + 1i*(w+wa)));

F = @(w) abs(Fu(w));

w = linspace(0, 2, 10000);
Fw = F(w);

[~, kwm] = max(Fw);
wm = w(kwm);

Fw1 = Fw(1:kwm);
kw1 = length(Fw1(Fw1 <= max(Fw)/sqrt(2)));
w1 = w(kw1);

Fw2 = Fw(kwm:end);
kw2 = length(Fw) - length(Fw2(Fw2 <= max(Fw)/sqrt(2)));
w2 = w(kw2);


fig = figure;
ax = axes(fig);
hold(ax, 'on');
plot(ax, w, Fw);
xlabel(ax, '\omega', 'FontSize', 15);
ylabel(ax, '|F[u](\omega)|', 'FontSize', 12);
ylim([0, 11]);

plot(ax, [w1, w2], [max(Fw)/sqrt(2), max(Fw)/sqrt(2)], 'b-*');
plot(ax, wm, max(Fw), 'r*');

text(ax, (w1+w2)/2 - 0.07, max(Fw)/sqrt(2) - 0.4, '\Delta\omega', 'FontSize', 15);
text(ax, wm + 0.07, max(Fw) + 0, '\omega_{max}', 'FontSize', 15);
text(ax, 0.1, 10, {'\omega_0 = 1', '\zeta = 5%'}, 'FontSize', 12);