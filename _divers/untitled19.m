phi = exp(1i*pi*[0 1 0.9 0.1]).';

Theta = [];
Theta0 = linspace(0, 2*pi, 5000);
for theta0 = Theta0
    phi = phi / exp(1i*theta0);
    phi = phi*exp(1i*pi/4);
    phir = real(phi);
    phii = imag(phi);
    
    theta = atan2((phir.'*phii), (phir.'*phir));
    
    Theta(end+1) = theta ;
end

figure;
plot(Theta0/pi, (Theta - Theta0)/pi);