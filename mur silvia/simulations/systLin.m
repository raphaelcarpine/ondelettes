classdef systLin
    %SYSTLIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M {mustBeNumeric}
        K {mustBeNumeric}
        C {mustBeNumeric}
        n {mustBeNumeric}
    end
    
    methods
        function obj = systLin(M, K, C)
            obj.M = M;
            obj.K = K;
            if nargin < 3
                C = zeros(size(M));
            end
            obj.C = C;
            obj.n = size(M, 1);
            
            if size(M, 1) ~= size(M, 2) || ~all(size(M) == size(K)) || ~all(size(M) == size(C))
                error('');
            end
            for Mat = {M, K, C}
                if ~all(transpose(Mat{1}) == Mat{1})
                    error('');
                end
            end
        end
        
        function [poles, shapes] = complexModes(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Mat = [zeros(obj.n), eye(obj.n);
                -obj.M\obj.K, -obj.M\obj.C];
            [shapes0, poles0] = eig(Mat);
            poles0 = diag(poles0);
            shapes0 = shapes0(1:obj.n, :);
            
            poles = [];
            shapes = [];
            for k = 1:length(poles0)
                if imag(poles0(k)) >= 0
                    poles(end+1) = poles0(k);
                    shapes = [shapes, normalisation(obj, shapes0(:, k))];
                end
            end
            
            % tri croissant
            [~, I] = sort(poles/1i);
            poles = poles(I);
            shapes = shapes(:, I);
        end
        
        function [poles, shapes] = normalModes(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [poles, shapes] = complexModes(systLin(obj.M, obj.K));
            shapes = real(shapes);
        end
        
        function [Mb, Kb, Cb] = modalDamping(obj)
            [~, shapes] = normalModes(obj);
            Matb = shapes * (shapes.'*obj.M*shapes)^(-1/2);
            Mb = Matb' * obj.M * Matb;
            Kb = Matb' * obj.K * Matb;
            Cb = Matb' * obj.C * Matb;
        end
        
        function x = response(obj, f, dt, ddl)
            nT = size(f, 2);
            if nargin == 4
                f0 = zeros(obj.n, nT);
                f0(ddl, :) = f;
                f = f0;
            end
            
            % vecteur i\omega
            iw = (0:nT-1) * (2i*pi/(nT*dt));
            for k = 1:floor((nT-1)/2)
                iw(end+1-k) = -k * (2i*pi/(nT*dt));
            end
            
            FFTf = fft(f, nT, 2);
            FFTx = nan(size(FFTf));
            for k = 1:nT % * fct transfert
                FFTx(:, k) = (iw(k)^2*obj.M + iw(k)*obj.C + obj.K)\FFTf(:, k);
            end
            
            x = ifft(FFTx, nT, 2);
        end
    end
    
    methods (Access = private)
        function shape = normalisation(obj, shape)
            shape = shape / sqrt(transpose(shape) * shape);
            [~, k] = max(abs(shape));
            shape = shape * sign(real(shape(k)));
        end
    end
end

