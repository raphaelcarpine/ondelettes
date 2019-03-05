classdef TMD
    %TMD 
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        m0 % masse
        m1 % fonction de masse équivalente : m1(x)x"+m0x0"=f(x,x')
        f % fonction de forces : m1(x)x"+m0x0"=f(x,x')
        f0 % force appliquée sur le bati : F0 = f0{1}(x,x')x0" + f0{2}(x,x') 
    end
    
    methods
        function obj = TMD(m0, m1, f, f0)
            %TMD Construct an instance of this class
            %   Detailed explanation goes here
            obj.m0 = m0;
            obj.m1 = m1;
            obj.f = f;
            obj.f0 = f0;
        end
        
        function [t, x] = reponseLibre(obj, x0, v0, T)
            %reponseLibre Summary of this method goes here
            %   Detailed explanation goes here
            F = @(t, X) [X(2); 1/obj.m1(X(1))*obj.f(X(1), X(2))];
            [t, x] = differentialEq([x0 v0], F, T, true);
            animate(obj, t, x);
        end
        
        function amortissementFreq(obj, ampl, freqs, periodesRegimeTransitoire)
            %reponseLibre Summary of this method goes here
            %   Detailed explanation goes here
            puissances = zeros(1, length(freqs));
            wait = waitbar(0, sprintf('0/%d', length(freqs)), 'Name', "Calcul de l'amortissement");
            for i = 1:length(freqs)
                freq = freqs(i);
                omega = 2*pi*freq;
                T = periodesRegimeTransitoire / freq;
                F = @(t, X) [X(2); 1/obj.m1(X(1))*obj.f(X(1), X(2)) + obj.m0/obj.m1(X(1))*omega^2*ampl*sin(omega*t)];
                [~, X] = differentialEq([0 0], F, T, false);
                T = 1/freq;
                [t, X] = differentialEq(X(end,:), F, T, false);
                puissance = omega*ampl*cos(omega*t) .* obj.f(X(:,1), X(:,2));
%                 figure;
%                 plot(t, puissance);
                I = 0;
                for j = 1:length(t)-1
                    I = I + (puissance(j)+puissance(j+1))/2 * (t(j+1)-t(j));
                end
                puissances(i) = freq * I;
                waitbar(i/length(freqs), wait, sprintf('%d/%d',i,length(freqs)));
            end
            close(wait);
            figure();
            plot(freqs, puissances);
            xlabel('f');
            ylabel('P');
        end
    end
end

