w1 = 1;
w2 = 1.001;

M = eye(2);
K = diag([w1^2, w2^2]);
C = 0.01*[1, -1; -1, 1];


syst = systLin(M, K, C);


[poles, shapes] = syst.normalModes()

[poles, shapes] = syst.complexModes()

phi1 = shapes(:, 1);
phi2 = shapes(:, 2);

figure;
polarplot([0, phi1(1)]);
hold on
polarplot([0, phi1(2)]);

figure;
polarplot([0, phi2(1)]);
hold on
polarplot([0, phi2(2)]);

return

%%

for theta = linspace(pi, 0, 20)
    phi = [exp(1i*theta); exp(-1i*theta)];
    phi = phi / sqrt(phi.'*phi);
    figure;
    polarplot([0, phi(1)]);
    hold on
    polarplot([0, phi(2)]);
end