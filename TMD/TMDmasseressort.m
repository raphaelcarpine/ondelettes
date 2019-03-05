classdef TMDmasseressort < TMD
    %TMDmasseressort Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        m % masse
        k % constante de raideur
        a % fonction d'amortissement a(x, x')
    end
    
    methods
        function obj = TMDmasseressort(m, k, a)
            %TMDmasseressort Construct an instance of this class
            %   Detailed explanation goes here
            obj@TMD(m, @(x) m, @(x,v) -k*x-a(x,v), {@(x,v) 0, @(x,v) k*x+a(x,v)});
            obj.m = m;
            obj.k = k;
            obj.a = a;
        end
    end
end

