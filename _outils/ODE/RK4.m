function y = RK4(f, t, y0)
%RK4 Runge-Kutta 4
%   y' = f(t, y)
%   y column

[initWaitBar, updateWaitBar, closeWaitBar] = getWaitBar(length(t)-1,...
    'displayTime', 0, 'windowTitle', 'Integrating ODE (RK4)');


y = nan(length(y0), length(t));
y(:, 1) = y0;

initWaitBar();
for kt = 1:length(t)-1
    tk = t(kt);
    yk = y(:, kt);
    h = t(kt+1) - tk;
    
    k1 = f(tk, yk);
    k2 = f(tk + h/2, yk + h/2*k1);
    k3 = f(tk + h/2, yk + h/2*k2);
    k4 = f(tk + h, yk + h*k3);
    
    y(:, kt+1) = yk + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    updateWaitBar(kt+1);
end
closeWaitBar();


end

