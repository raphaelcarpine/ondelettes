classdef LongInt
    %LONGINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant, Access = private)
        kmax = int64(3037000499); % 0 <= K(i) < kmax
        logkmax = log(3037000499);
    end
    
    properties (Constant)
        Kmax = LongInt(LongInt.kmax);
    end
    
    properties
        K % int array, int64
        s % sign, bool
    end
    
    methods
        function obj = LongInt(k)
            %LONGINT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 0
                obj.s = k >= 0;
                obj.K = [mod(abs(int64(k)), LongInt.kmax), idivide(abs(int64(k)), LongInt.kmax)];
                if obj.K(2) >= LongInt.kmax
                    obj.K(2) = obj.K(2) - LongInt.kmax;
                    obj.K(3) = 1;
                end
                while ~isempty(obj.K) &&  obj.K(end) == 0
                    obj.K(end) = [];
                end
            end
        end
        
        function disp(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            n = length(obj.K);
            Ksym = sym(obj.K);
            ksym = sym(LongInt.kmax);
            disp(sum(Ksym .* ksym.^(0:n-1)));
            fprintf('(%d int64 digits)\n', length(obj.K));
        end
        
        function obj2 = copy(obj)
            obj2 = LongInt;
            obj2.s = obj.s;
            obj2.K = obj.K;
        end
        
        function obj = trunc(obj, n)
            if n > 0
                if n < length(obj.K)
                    obj.K(1:n) = [];
                else
                    obj.K = int64([]);
                end
            elseif n < 0
                obj.K = [int64(zeros(1, -n)), obj.K];
            end
        end
        
        function obj3 = plus(obj1, obj2)
            if obj1.s == obj2.s
                n1 = length(obj1.K);
                n2 = length(obj2.K);
                n = max(n1, n2) + 1;
                K3 = [obj1.K, zeros(1, n-n1)] + [obj2.K, zeros(1, n-n2)];
                for i = 1:n-1
                    if K3(i) >= LongInt.kmax
                        K3(i) = K3(i) - LongInt.kmax;
                        K3(i+1) = K3(i+1) + 1;
                    end
                end
                if K3(end) == 0
                    K3(end) = [];
                end
                obj3 = LongInt;
                obj3.s = obj1.s;
                obj3.K = K3;
            else
                obj2.s = ~obj2.s;
                obj3 = obj1 - obj2;
                obj2.s = ~obj2.s;
            end
        end
        
        function obj3 = minus(obj1, obj2)
            if obj1.s == obj2.s
                obj3 = LongInt;
                obj3.s = obj1.s;
                if xor(obj1 > obj2, obj1.s) % if abs(obj1) <= abs(obj2)
                    [obj2, obj1] = deal(obj1, obj2);
                    obj3.s = ~obj3.s;
                end
                n1 = length(obj1.K);
                n2 = length(obj2.K);
                n = n1;
                K3 = obj1.K - [obj2.K, zeros(1, n-n2)];
                for i = 1:n-1
                    if K3(i) < 0
                        K3(i) = K3(i) + LongInt.kmax;
                        K3(i+1) = K3(i+1) - 1;
                    end
                end
                while ~isempty(K3) && K3(end) == 0
                    K3(end) = [];
                end
                obj3.K = K3;
            else
                obj2.s = ~obj2.s;
                obj3 = obj1 + obj2;
                obj2.s = ~obj2.s;
            end
        end
        
        function obj3 = times(obj1, obj2)
            if ~isa(obj1, 'LongInt')
                [obj2, obj1] = deal(obj1, obj2);
            end
            if isa(obj2, 'LongInt')
                obj3 = LongInt;
                obj3.s = ~xor(obj1.s, obj2.s);
                n1 = length(obj1.K);
                n2 = length(obj2.K);
                n = n1 + n2;
                K3 = int64(zeros(1, n));
                for i = 1:n1
                    for j = 1:n2
                        M = obj1.K(i) * obj2.K(j);
                        r = mod(M, LongInt.kmax);
                        d = (M-r) ./ LongInt.kmax;
                        K3(i+j-1) = K3(i+j-1) + r;
                        n0 = i+j-1;
                        while K3(n0) >= LongInt.kmax
                            K3(n0) = K3(n0) - LongInt.kmax;
                            K3(n0+1) = K3(n0+1) + 1;
                            n0 = n0 + 1;
                        end
                        K3(i+j) = K3(i+j) + d;
                        n0 = i+j;
                        while K3(n0) >= LongInt.kmax
                            K3(n0) = K3(n0) - LongInt.kmax;
                            K3(n0+1) = K3(n0+1) + 1;
                            n0 = n0 + 1;
                        end
                    end
                end
                while ~isempty(K3) &&  K3(end) == 0
                    K3(end) = [];
                end
                obj3.K = K3;
            else
                x = double(obj2);
                obj3 = obj1.copy();
                if x < 0
                    obj3.s = ~obj3.s;
                    x = abs(x);
                end
                nx = floor(x);
                rx = mod(x, 1);
                r = 0;
                for i = length(obj3.K):-1:1
                    dk = double(obj3.K(i))*rx + double(LongInt.kmax)*r;
                    r = mod(dk, 1);
                    if i > 1
                        obj3.K(i) = int64(floor(dk));
                    else
                        obj3.K(i) = int64(round(dk));
                    end
                end
                while ~isempty(obj3.K) && obj3.K(end) == 0
                    obj3.K(end) = [];
                end
                obj3 = obj3 + obj1 .* LongInt(nx);
            end
        end
        
        function c = gt(obj1, obj2)
            if obj1.s ~= obj2.s
                c = obj1.s;
            else
                if length(obj1.K) ~= length(obj2.K)
                    c = length(obj1.K) > length(obj2.K);
                else
                    c = false;
                    for i = length(obj1.K):-1:1
                        if obj1.K(i) > obj2.K(i)
                            c = true;
                            break
                        elseif obj1.K(i) < obj2.K(i)
                            c = false;
                            break
                        end
                    end
                end
                c = xor(c, ~obj1.s);
            end
        end
        
        function x = double(obj)
            x = sum(exp(log(double(obj.K)) + LongInt.logkmax*(0:length(obj.K)-1)));
        end
        
        function x = divideDouble(obj1, obj2)
            K1 = flip(obj1.K);
            K2 = flip(obj2.K);
            x = sum(exp(log(double(K1)) - LongInt.logkmax*(0:length(K1)-1)));
            x = x / sum(exp(log(double(K2)) - LongInt.logkmax*(0:length(K2)-1)));
            x = x * exp((length(obj1.K) - length(obj2.K)) * LongInt.logkmax);
        end
    end
    
    methods (Static)
        function Z = zeros(n1, n2)
            Z = LongInt(0);
            for i = 1:n1
                for j = 1:n2
                    Z(i, j) = LongInt(0);
                end
            end
        end
    end
end

