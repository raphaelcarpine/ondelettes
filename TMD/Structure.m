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
        
        function [t, X] = reponseForceeSansTMD(obj, x0, v0, T, ext, visible)
            %reponseLibreSansTMD forçage sur le premier ddl seulement
            %   Detailed explanation goes here
            switch nargin
                case 5
                    visible = false;
            end
            
            function Y = F(t, X)
                x = X(1:obj.ddl);
                v = X(obj.ddl+1:2*obj.ddl);
                xext = zeros(length(x), 1);
                vext = zeros(length(x), 1);
                Ext = ext(t);
                xext(1) = Ext(1);
                vext(1) = Ext(2);
                
                forces = -obj.K*(x-xext) - obj.A(x-xext, v-vext);
                Y = zeros(2*obj.ddl, 1);
                Y(1:obj.ddl) = v;
                Y(obj.ddl+1:2*obj.ddl) = obj.M\forces;
            end
            [t, X] = differentialEq([x0 v0].', @F, T, visible);
            if visible
                animate(obj, t, X);
            end
        end
        
        function [t, X] = reponseForceeAvecTMD(obj, x0, v0, T, ext, visible)
            %reponseLibreSansTMD Summary of this method goes here
            %   Detailed explanation goes here
            switch nargin
                case 5
                    visible = false;
            end
            
            nddl = obj.ddl;
            ntmd = length(obj.TMDs);
            function Y = F(t, X)
                x = X(1:nddl);
                xtmd = X(nddl+1:nddl+ntmd);
                v = X(nddl+ntmd+1:2*nddl+ntmd);
                vtmd = X(2*nddl+ntmd+1:2*nddl+2*ntmd);
                xext = zeros(length(x), 1);
                vext = zeros(length(x), 1);
                Ext = ext(t);
                xext(1) = Ext(1);
                vext(1) = Ext(2);
                
                forces = -obj.K*(x-xext) - obj.A(x-xext, v-vext);
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
            [t, X] = differentialEq(X0, @F, T, visible);
            if visible
                animate(obj, t, X);
            end
        end
        
        function [t, X] = reponseLibre(obj, x0, v0, T, avecTMD)
            %reponseLibreSansTMD Summary of this method goes here
            %   Detailed explanation goes here
            if avecTMD
                reponseForcee = @(varargin) obj.reponseForceeAvecTMD(varargin{:});
            else
                reponseForcee = @(varargin) obj.reponseForceeSansTMD(varargin{:});
            end
            
            [t, X] = reponseForcee(x0, v0, T, @(x, v) [0 0], true);
        end
        
        function diagrammeBode(obj, ddl, ampl, freqs, tempsRegimeTransitoire, avecTMD)
            %reponseLibre Summary of this method goes here
            %   Detailed explanation goes here
            if avecTMD
                reponseForcee = @(varargin) obj.reponseForceeAvecTMD(varargin{:});
            else
                reponseForcee = @(varargin) obj.reponseForceeSansTMD(varargin{:});
            end
            
            amplitudes = zeros(1, length(freqs));
            wait = waitbar(0, sprintf('0/%d', length(freqs)), 'Name', "Calcul du diagramme de Bode");
            for i = 1:length(freqs)
                freq = freqs(i);
                omega = 2*pi*freq;
                T = tempsRegimeTransitoire + 1/freq;
                x0 = zeros(1, obj.ddl);
                v0 = zeros(1, obj.ddl);
                [t, X] = reponseForcee(x0, v0, T, @(t) [ampl*sin(omega*t), ampl*omega*cos(omega*t)]);
                
                n = length(t);
                while t(n) >= tempsRegimeTransitoire
                    n = n-1;
                end
                
                amplitudes(i) = (max(X(n+1:end, ddl))-min(X(n+1:end, ddl)))/2;
                waitbar(i/length(freqs), wait, sprintf('%d/%d',i,length(freqs)));
            end
            close(wait);
            figure;
            loglog(freqs, amplitudes/ampl);
            set(gca, 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':');
            xlabel('f');
            ylabel('H');
        end
    end
end

