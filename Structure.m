classdef Structure
    %Structure Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        ddl % nombre de ddl
        M % matrice de masse, DOIT ETRE DIAGONALE
        K % matrice de raideur
        A % fonction d'amortissement A(x,x')
    end
    
    properties (SetAccess = private)
        TMDs % cellule contenant les TMDs de la structure, avec leur ddl {{TMD1, ddl1},...
    end
    
    methods
        function obj = Structure(M, K, A, TMDs)
            %Structure Construct an instance of this class
            %   Detailed explanation goes here
            obj.M = M;
            obj.K = K;
            obj.A = A;
            obj.ddl = length(M);
            obj.TMDs = TMDs;
        end
        
        function [t, x] = reponseLibreSansTMD(obj, x0, v0, T)
            %reponseLibreSansTMD Summary of this method goes here
            %   Detailed explanation goes here
            function Y = F(~, X)
                x = X(1:obj.ddl);
                v = X(obj.ddl+1:2*obj.ddl);
                forces = -obj.K*x - obj.A(x, v);
                Y = zeros(2*obj.ddl, 1);
                Y(1:obj.ddl) = v;
                Y(obj.ddl+1:2*obj.ddl) = obj.M\forces;
            end
            [t, x] = differentialEq([x0 v0].', @F, T, true);
        end
        
        function [t, x] = reponseLibreAvecTMD(obj, x0, v0, T)
            %reponseLibreSansTMD Summary of this method goes here
            %   Detailed explanation goes here
            nddl = obj.ddl;
            ntmd = length(obj.TMDs);
            function Y = F(~, X)
                x = X(1:nddl);
                xtmd = X(nddl+1:nddl+ntmd);
                v = X(nddl+ntmd+1:2*nddl+ntmd);
                vtmd = X(2*nddl+ntmd+1:2*nddl+2*ntmd);
                
                forces = -obj.K*x - obj.A(x, v);
                forcestmd =  zeros(ntmd,1);
                for itmd=1:ntmd
                    tmd = obj.TMDs{itmd}{1};
                    ddltmd = obj.TMDs{itmd}{2};
                    forces(ddltmd) = forces(ddltmd) + tmd.f0{2}(xtmd(itmd), vtmd(itmd));
                    forcestmd(itmd) = tmd.f(xtmd(itmd), vtmd(itmd));
                end
                
                M2 = zeros(nddl+ntmd);
                M2(1:nddl, 1:nddl) = obj.M;
                for itmd=1:ntmd
                    tmd = obj.TMDs{itmd}{1};
                    ddltmd = obj.TMDs{itmd}{2};
                    M2(ddltmd, ddltmd) = M2(ddltmd, ddltmd) - tmd.f0{1}(xtmd(itmd), vtmd(itmd));
                    M2(nddl+itmd, ddltmd) = tmd.m0;
                    M2(nddl+itmd, nddl+itmd) = tmd.m1(xtmd(itmd));
                end
                
                Y = zeros(2*nddl+2*ntmd, 1);
                Y(1:nddl+ntmd) = [v;vtmd];
                Y(nddl+ntmd+1:2*nddl+2*ntmd) = M2\[forces;forcestmd];
            end
            
            v0tmd = zeros(1,ntmd);
            for jtmd=1:ntmd
                tmd2 = obj.TMDs{jtmd}{1};
                ddltmd2 = obj.TMDs{jtmd}{2};
                v0tmd(jtmd) = -v0(ddltmd2)*tmd2.m0/tmd2.m1(0);
            end
            X0 = [x0 zeros(1,ntmd) v0 v0tmd].';
            [t, x] = differentialEq(X0, @F, T, true);
        end
    end
end

