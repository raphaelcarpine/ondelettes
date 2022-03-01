function ridge = SingleRidgeExtract(T, freqs, CWT, MotherWavelet, Q, ct, ridgeContinuity, slopeTimeConst)
%SINGLERIDGEEXTRACT Summary of this function goes here
%   Detailed explanation goes here

dt = (T(end)-T(1))/(length(T)-1);

% ridge computation
switch ridgeContinuity
    case 'none'
        freqs = [nan, freqs, nan];
        CWT = [zeros(1, size(CWT, 2)); CWT; zeros(1, size(CWT, 2))];
        [~, Fridge] = max(abs(CWT), [], 1);
        [Fridge, Aridge] = localMax3Points(freqs([Fridge-1; Fridge; Fridge+1]),...
            CWT([Fridge-1; Fridge; Fridge+1] + [1;1;1] * (0:size(CWT, 2)-1)*size(CWT, 1)));
    case 'simple'
        freqs = [nan, freqs, nan];
        CWT = [zeros(1, size(CWT, 2)); CWT; zeros(1, size(CWT, 2))];
        Fridge = nan(1, size(CWT, 2));
        [~, Fridge(1)] = max(abs(CWT(:, 1)));
        localMax = abs(CWT(1:end-1, :)) < abs(CWT(2:end, :));
        localMax = localMax(1:end-1, :) & ~localMax(2:end, :);
        for kt = 2:size(CWT, 2)
            localMaxFreq = find(localMax(:, kt)) + 1;
            [~, closestLocalMax] = min(abs(localMaxFreq - Fridge(kt-1)));
            Fridge(kt) = localMaxFreq(closestLocalMax);
        end
        [Fridge, Aridge] = localMax3Points(freqs([Fridge-1; Fridge; Fridge+1]),...
            CWT([Fridge-1; Fridge; Fridge+1] + [1;1;1] * (0:size(CWT, 2)-1)*size(CWT, 1)));
    case 'reverse'
        freqs = [nan, freqs, nan];
        CWT = [zeros(1, size(CWT, 2)); CWT; zeros(1, size(CWT, 2))];
        Fridge = nan(1, size(CWT, 2));
        [~, Fridge(end)] = max(abs(CWT(:, end)));
        localMax = abs(CWT(1:end-1, :)) < abs(CWT(2:end, :));
        localMax = localMax(1:end-1, :) & ~localMax(2:end, :);
        for kt = size(CWT, 2)-1:-1:1
            localMaxFreq = find(localMax(:, kt)) + 1;
            [~, closestLocalMax] = min(abs(localMaxFreq - Fridge(kt+1)));
            Fridge(kt) = localMaxFreq(closestLocalMax);
        end
        [Fridge, Aridge] = localMax3Points(freqs([Fridge-1; Fridge; Fridge+1]),...
            CWT([Fridge-1; Fridge; Fridge+1] + [1;1;1] * (0:size(CWT, 2)-1)*size(CWT, 1)));
    case 'slope' % slope penalization, 10.1109/78.640725 Eq.(4)
        % frequency increment
        df = (freqs(end) - freqs(1))/(length(freqs)-1);
        if any(abs(diff(freqs)/df - 1) > 1e-4)
            error('frequency increment must be constant');
        end
        
        % signal mean square
        [~, Fridge0] = max(abs(CWT), [], 1);
        Aridge0 = CWT(Fridge0 + (0:size(CWT, 2)-1)*size(CWT, 1));
        A2ref = mean(abs(Aridge0).^2);
        Fridge0(5:end) = freqs(Fridge0(5:end));
        
        % first ridge values correction
        [Fridge0(1:4), Aridge0(1:4)] = localMax3Points(freqs([Fridge0(1:4)-1; Fridge0(1:4); Fridge0(1:4)+1]),...
            CWT([Fridge0(1:4)-1; Fridge0(1:4); Fridge0(1:4)+1] + [1;1;1] * (0:3)*size(CWT, 1)));
        
        % ODE
        lambda = A2ref * slopeTimeConst^2;
        D = @(t, phi) [phi(2); - exp(phi(1))/(2*lambda)...
            * doubleQuadInterpDeriv(round((t-T(1))/dt)+1, exp(phi(1)))];
        phi0 = [log(Fridge0(1)); (log(Fridge0(2))-log(Fridge0(1)))/dt];
        phi = RK4(D, T, phi0);
        phi = phi(1, :);
        
        % ridge freq and ampl
        Fridge = exp(phi);
        Aridge = nan(size(Fridge));
        for kt = 1:length(Fridge)
            Aridge(kt) = doubleQuadInterpAmpl(kt, Fridge(kt));
        end
    case 'slope2' % slope penalization, 10.1109/78.640725 Eq.(4)
        % frequency increment
        df = (freqs(end) - freqs(1))/(length(freqs)-1);
        if any(abs(diff(freqs)/df - 1) > 1e-4)
            error('frequency increment must be constant');
        end
        
        % signal mean square
        [~, Fridge0] = max(abs(CWT), [], 1);
        Aridge0 = CWT(Fridge0 + (0:size(CWT, 2)-1)*size(CWT, 1));
        A2ref = mean(abs(Aridge0).^2);
%         Fridge0(5:end) = freqs(Fridge0(5:end));
        
        % first ridge values correction
%         [Fridge0(1:4), Aridge0(1:4)] = localMax3Points(freqs([Fridge0(1:4)-1; Fridge0(1:4); Fridge0(1:4)+1]),...
%             CWT([Fridge0(1:4)-1; Fridge0(1:4); Fridge0(1:4)+1] + [1;1;1] * (0:3)*size(CWT, 1)));
        [Fridge0, Aridge0] = localMax3Points(freqs([Fridge0-1; Fridge0; Fridge0+1]),...
            CWT([Fridge0-1; Fridge0; Fridge0+1] + [1;1;1] * (0:size(CWT, 2)-1)*size(CWT, 1)));
        
        % ODE
        lambda = A2ref * slopeTimeConst^2;
        D = @(t, phi) [phi(2); - exp(phi(1))/(2*lambda)...
            * doubleQuadInterpDeriv(round((t-T(1))/dt)+1, exp(phi(1)))];
        phi0 = log(Fridge0(1));
        phi1 = log(Fridge0(end));
        dphi0u = abs((log(Fridge0(2))-log(Fridge0(1)))/dt);
        dphi0l = -abs((log(Fridge0(2))-log(Fridge0(1)))/dt);
        while true
            phi = RK4(D, T, [phi0; dphi0u]);
            if phi(end) >= phi1
                break
            else
               dphi0u = 2*dphi0u;
            end
        end
        while true
            phi = RK4(D, T, [phi0; dphi0l]);
            if phi(end) <= phi1
                break
            else
               dphi0l = 2*dphi0l;
            end
        end
        for k0 = 1:1000
            dphi0 = (dphi0l + dphi0u)/2;
            phi = RK4(D, T, [phi0; dphi0]);
            disp(phi(1, end)-phi1);
            if phi(1, end) >= phi1
                dphi0u = dphi0;
            else
                dphi0l = dphi0;
            end
        end
        phi = phi(1, :);
        
        % ridge freq and ampl
        Fridge = exp(phi);
        Aridge = nan(size(Fridge));
        for kt = 1:length(Fridge)
            Aridge(kt) = doubleQuadInterpAmpl(kt, Fridge(kt));
        end
    case 'slope3' % slope penalization, 10.1109/78.640725 Eq.(4)
        % frequency increment
        df = (freqs(end) - freqs(1))/(length(freqs)-1);
        if any(abs(diff(freqs)/df - 1) > 1e-4)
            error('frequency increment must be constant');
        end
        
        % signal mean square
        [~, Fridge0] = max(abs(CWT), [], 1);
        Aridge0 = CWT(Fridge0 + (0:size(CWT, 2)-1)*size(CWT, 1));
        A2ref = mean(abs(Aridge0).^2);
        [Fridge0, ~] = localMax3Points(freqs([Fridge0-1; Fridge0; Fridge0+1]),...
            CWT([Fridge0-1; Fridge0; Fridge0+1] + [1;1;1] * (0:size(CWT, 2)-1)*size(CWT, 1)));
        
        % ODE
        lambda = A2ref * slopeTimeConst^2;
        phi = log(Fridge0(2:end-1)).';
        DD = zeros(length(phi));
        DD(1, 1:2) = [-2 1];
        DD(end, end-1:end) = [1 -2];
        for kl = 2:length(phi)-1
            DD(kl, kl-1:kl+1) = [1 -2 1];
        end
        figure;
        for k = 1:100
            plot(T(2:end-1), exp(phi));
            disp(k);
            drawnow
            input(' ');
            for kt = 1:length(phi)
                phi(kt) = - dt^2 * exp(phi(kt))/(2*lambda) *...
                    doubleQuadInterpDeriv(round((T(kt)-T(1))/dt)+1, exp(phi(kt)));
            end
            phi([1, end]) = phi([1, end]) - log(Fridge0([1, end])).';
            phi = DD \ phi;
        end
        
        % ridge freq and ampl
        Fridge = [Fridge0(1), exp(phi), Fridge0(end)];
        Aridge = nan(size(Fridge));
        for kt = 1:length(Fridge)
            Aridge(kt) = doubleQuadInterpAmpl(kt, Fridge(kt));
        end
    case 'curv' % slope penalization, 10.1109/78.640725 Eq.(4)
        % frequency increment
        df = (freqs(end) - freqs(1))/(length(freqs)-1);
        if any(abs(diff(freqs)/df - 1) > 1e-4)
            error('frequency increment must be constant');
        end
        
        % signal mean square
        [~, Fridge0] = max(abs(CWT), [], 1);
        Aridge0 = CWT(Fridge0 + (0:size(CWT, 2)-1)*size(CWT, 1));
        A2ref = mean(abs(Aridge0).^2);
        Fridge0(5:end) = freqs(Fridge0(5:end));
        
        % first ridge values correction
        [Fridge0(1:4), Aridge0(1:4)] = localMax3Points(freqs([Fridge0(1:4)-1; Fridge0(1:4); Fridge0(1:4)+1]),...
            CWT([Fridge0(1:4)-1; Fridge0(1:4); Fridge0(1:4)+1] + [1;1;1] * (0:3)*size(CWT, 1)));
        
        % ODE
        mu = A2ref * slopeTimeConst^4;
        D = @(t, phi) [phi(2:4); exp(phi(1))/(2*mu)...
            * doubleQuadInterpDeriv(round((t-T(1))/dt)+1, exp(phi(1)))...
            + phi(3)/slopeTimeConst^2*10];
        phi0 = [log(Fridge0(1)); (log(Fridge0(2))-log(Fridge0(1)))/dt;...
            (log(Fridge0(3))-2*log(Fridge0(2))+log(Fridge0(1)))/dt^2; 0];
        phi0 = [log(10); 0; 0; 0];
        phi = RK4(D, T, phi0);
        phi = phi(1, :);
        
        % ridge freq and ampl
        Fridge = exp(phi);
        Aridge = nan(size(Fridge));
        for kt = 1:length(Fridge)
            Aridge(kt) = doubleQuadInterpAmpl(kt, Fridge(kt));
        end
    case 'curv2' % slope penalization, 10.1109/78.640725 Eq.(4)
        % frequency increment
        df = (freqs(end) - freqs(1))/(length(freqs)-1);
        if any(abs(diff(freqs)/df - 1) > 1e-4)
            error('frequency increment must be constant');
        end
        
        % signal mean square
        [~, Fridge0] = max(abs(CWT), [], 1);
        Aridge0 = CWT(Fridge0 + (0:size(CWT, 2)-1)*size(CWT, 1));
        A2ref = mean(abs(Aridge0).^2);
        Fridge0(5:end) = freqs(Fridge0(5:end));
        
        % first ridge values correction
        [Fridge0(1:4), Aridge0(1:4)] = localMax3Points(freqs([Fridge0(1:4)-1; Fridge0(1:4); Fridge0(1:4)+1]),...
            CWT([Fridge0(1:4)-1; Fridge0(1:4); Fridge0(1:4)+1] + [1;1;1] * (0:3)*size(CWT, 1)));
        
        % ODE
        mu = A2ref * slopeTimeConst^4;
        D = @(t, phi) [phi(2:4); 1/(2*mu)...
            * doubleQuadInterpDeriv(round((t-T(1))/dt)+1, phi(1))...
            + phi(3)/slopeTimeConst^2/inf];
        phi0 = [Fridge0(1); (Fridge0(2)-Fridge0(1))/dt;...
            (Fridge0(3)-2*Fridge0(2)+Fridge0(1))/dt^2; 0];
        phi0 = [Fridge0(1); 0; 0; 0];
        phi = RK4(D, T, phi0);
        phi = phi(1, :);
        
        % ridge freq and ampl
        Fridge = phi;
        Aridge = nan(size(Fridge));
        for kt = 1:length(Fridge)
            Aridge(kt) = doubleQuadInterpAmpl(kt, Fridge(kt));
        end
    otherwise
        error(' ');
end

ridge.time = T;
ridge.val = Aridge;
ridge.freq = Fridge;

% edge effects
[~, DeltaT] = FTpsi_DeltaT(Q, MotherWavelet);
Ntedge = ceil(ct*DeltaT(freqs(1))/dt);
nt1 = 1;
for kt = min(Ntedge + 1, length(T)) : -1 : 1
    if ct*DeltaT(ridge.freq(kt)) > (kt-1) * dt
        nt1 = kt;
        break
    end
end
nt2 = length(T);
for kt = min(Ntedge, length(T)-1) : -1 : 0
    if ct*DeltaT(ridge.freq(length(T) - kt)) > kt * dt
        nt2 = length(T)- kt;
        break
    end
end

ridge.time = ridge.time(nt1:nt2);
ridge.val = ridge.val(nt1:nt2);
ridge.freq = ridge.freq(nt1:nt2);


%%

    function dT2 = doubleQuadInterpDeriv(kt, f)
        % frequencies indexes
        if f < freqs(1)
            dT2 = abs(CWT(1, kt))^2/df;
            return
        elseif f >= freqs(end)
            dT2 = - abs(CWT(end, kt))^2/df;
            return
        elseif f <= freqs(2)
            Kf1 = [1 2 3];
            Kf2 = Kf1;
        elseif f > freqs(end-1)
            Kf1 = length(freqs)-1 + [-1 0 1];
            Kf2 = Kf1;
        else
            Kf1 = floor((f - freqs(1))/df) + 1 + [-1 0 1];
            Kf2 = Kf1 + 1;
        end
        
        % intercepts
        f1 = freqs(Kf1);
        f2 = freqs(Kf2);
        T21 = CWT(Kf1, kt).';
        T22 = CWT(Kf2, kt).';
        T21 = abs(T21).^2; % square
        T22 = abs(T22).^2;
        
        % quadratic computation
        dT21 = sum(T21 .* ([1/2 -1 1/2]/df^2) .* (-sum(f1)*[1 1 1] + f1 + 2*f));
        dT22 = sum(T22 .* ([1/2 -1 1/2]/df^2) .* (-sum(f2)*[1 1 1] + f2 + 2*f));
        dT2 = (dT21 + dT22)/2;
    end

    function A = doubleQuadInterpAmpl(kt, f)
        % frequencies indexes
        if f < freqs(1) || f > freqs(end)
            A = nan;
            return
        elseif f < freqs(2)
            Kf1 = [1 2 3];
            Kf2 = Kf1;
        elseif f > freqs(end-1)
            Kf1 = length(freqs)-1 + [-1 0 1];
            Kf2 = Kf1;
        else
            Kf1 = floor((f - freqs(1))/df) + 1 + [-1 0 1];
            Kf2 = Kf1 + 1;
        end
        
        % intercepts
        f1 = freqs(Kf1);
        f2 = freqs(Kf2);
        T1 = CWT(Kf1, kt).';
        T2 = CWT(Kf2, kt).';
        
        % quadratic computation
        A1 = sum(T1 .* ([1/2 -1 1/2]/df^2) .* (prod(f-f1)*[1 1 1] ./ (f-f1)));
        A2 = sum(T2 .* ([1/2 -1 1/2]/df^2) .* (prod(f-f2)*[1 1 1] ./ (f-f2)));
        A = (A1 + A2)/2;
    end

end

