classdef TMDpendule < TMD
    %TMDpendule Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        g =  9.80665;
    end
    
    properties (SetAccess = immutable)
        I % moment d'inertie
        m % masse
        l % longueur
        a % fonction d'amortissement angulaire a(theta, theta')
    end
    
    methods
        function obj = TMDpendule(I, m, l, a)
            %TMDpendule Construct an instance of this class
            %   Detailed explanation goes here
            g = TMDpendule.g;
            m0 = m;
            m1 = @(theta) I/(l*cos(theta));
            f = @(theta, omega) -m*g*tan(theta) - a(theta, omega)/(l*cos(theta));
            f01 = @(theta, omega) (m*l*cos(theta))^2/I - m;
            f02 = @(theta, omega) (m*l)^2/I*g*cos(theta)*sin(theta) + m*l*cos(theta)/I*a(theta, omega) + m*l*omega^2*sin(theta);
            f0 = {f01, f02};
            
            obj@TMD(m0, m1, f, f0);
            obj.I = I;
            obj.m = m;
            obj.l = l;
            obj.a = a;
        end
    end
end

