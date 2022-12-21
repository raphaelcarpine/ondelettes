%%

Q = 5;

T = 40;
t = linspace(-T, T, 1000);
Psi = @(t) 1/(Q*sqrt(pi))*exp(1i*t-t.^2/(4*Q^2));
Psiab = @(a, b, t) Psi((t-b)/a);


fig = figure;
pltRe = plot(t, real(Psi(t)));
hold on
pltIm = plot(t, imag(Psi(t)));
xlabel('$\tau$', 'Interpreter', 'Latex', 'FontSize', 14);
ylabel('$\psi\left(\frac{\tau-b}{a}\right)$', 'Interpreter', 'Latex', 'FontSize', 14);
legend({'Re', 'Im'});
fig.Position(3:4) = [450 320];
xlim([-T T]);
% ylim([-0.15 0.15]);

%%

txtA = uicontrol('Style', 'text', 'String', 'a = 1.00', 'Position', [70 10 100 20], 'FontSize', 10);
txtB = uicontrol('Style', 'text', 'String', 'b = 0.0', 'Position', [220 10 100 20], 'FontSize', 10);

%% chemin

fps = 25;
C = [
    1    1    1    1      1    1    2    0.3  0.3  1;   % a
    0    0    10   -15    -10  -10  -10  -10  -10  0; % b
    0.5  0.5  1.25  0.25  0.5  0.5  1    0.5  0.5  0.5; % \delta t
    ];

updateFig(C(1, 1), C(2, 1), t, pltRe, pltIm, txtA, txtB, Psiab);
gif('test.gif', 'Resolution', 400, 'DelayTime', 1/fps, 'frame', gcf);

for k = 1:size(C, 2)
    kp1 = mod(k, size(C, 2)) + 1;
    nframes = round(C(3, k)*fps);
    coeff = linspace(0, 1, nframes+1);
    coeff(end) = [];
    A = (1-coeff)*C(1, k) + coeff*C(1, kp1);
    B = (1-coeff)*C(2, k) + coeff*C(2, kp1);
    for i = 1:length(coeff)
        updateFig(A(i), B(i), t, pltRe, pltIm, txtA, txtB, Psiab);
        drawnow;
        gif();
    end
end




function updateFig(a, b, t, pltRe, pltIm, txtA, txtB, Psiab)
    psi = Psiab(a, b, t);
    set(pltRe, 'YData', real(psi));
    set(pltIm, 'YData', imag(psi));
    set(txtA, 'String', sprintf('a = %.2f', a));
    set(txtB, 'String', sprintf('b = %.1f', b));
end