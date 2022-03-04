classdef LongInt
    %LONGINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        kmax = int64(3037000499); % 0 <= K(i) < kmax
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
                elseif obj.K(2) == 0
                    obj.K(2) = [];
                end
            end
        end
        
        function str = toString(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            n = length(obj.K);
            Ksym = sym(obj.K);
            ksym = sym(LongInt.kmax);
            str = char(sum(Ksym .* ksym.^(0:n-1)));
        end
        
        function obj2 = copy(obj)
            obj2 = LongInt;
            obj2.s = obj.s;
            obj2.K = obj.K;
        end
        
        function tronc(obj, n)
            obj.K(1:n) = [];
        end
    end
    
    methods (Static)
        function c = compare(obj1, obj2) % c=1 if obj1>obj2, c=-1 if obj1<obj2, c=0 if obj1=obj2
            if obj1.s ~= obj2.s
                c = 2*obj1.s - 1;
            else
                if length(obj1.K) ~= length(obj2.K)
                    c = 2*(length(obj1.K) > length(obj2.K)) - 1;
                else
                    c = 0;
                    for i = length(obj1.K):-1:1
                        if obj1.K(i) > obj2.K(i)
                            c = 1;
                            break
                        elseif obj1.K(i) < obj2.K(i)
                            c = -1;
                            break
                        end
                    end
                end
                c = c * (2*obj1.s - 1);
            end
        end
        
        function obj3 = add(obj1, obj2)
            if obj1.s == obj2.s
                n1 = length(obj1.K);
                n2 = length(obj2.K);
                n = max(n1, n2) + 1;
                K = [obj1.K, zeros(1, n-n1)] + [obj2.K, zeros(1, n-n2)];
                for i = 1:n-1
                    if K(i) >= LongInt.kmax
                        K(i) = K(i) - LongInt.kmax;
                        K(i+1) = K(i+1) + 1;
                    end
                end
                if K(end) == 0
                    K(end) = [];
                end
                obj3 = LongInt;
                obj3.s = obj1.s;
                obj3.K = K;
            else
                obj2.s = ~obj2.s;
                obj3 = LongInt.substract(obj1, obj2);
                obj2.s = ~obj2.s;
            end
        end
        
        function obj3 = substract(obj1, obj2)
            if obj1.s == obj2.s
                obj3 = LongInt;
                obj3.s = obj1.s;
                if LongInt.compare(obj1, obj2) * (2*obj1.s -1) < 0
                    [obj2, obj1] = deal(obj1, obj2);
                    obj3.s = ~obj3.s;
                end
                n1 = length(obj1.K);
                n2 = length(obj2.K);
                n = n1;
                K = obj1.K - [obj2.K, zeros(1, n-n2)];
                for i = 1:n-1
                    if K(i) < 0
                        K(i) = K(i) + LongInt.kmax;
                        K(i+1) = K(i+1) - 1;
                    end
                end
                if K(end) == 0
                    K(end) = [];
                end
                obj3.K = K;
            else
                obj2.s = ~obj2.s;
                obj3 = LongInt.add(obj1, obj2);
                obj2.s = ~obj2.s;
            end
        end
        
        function obj3 = multiply(obj1, obj2)
            obj3 = LongInt;
            obj3.s = ~xor(obj1.s, obj2.s);
            n1 = length(obj1.K);
            n2 = length(obj2.K);
            n = n1 + n2;
            K = int64(zeros(1, n));
            for i = 1:n1
                for j = 1:n2
                    M = obj1.K(i) * obj2.K(j);
                    r = mod(M, LongInt.kmax);
                    d = idivide(M, LongInt.kmax);
                    K(i+j-1) = K(i+j-1) + r;
                    n0 = i+j-1;
                    while K(n0) >= LongInt.kmax
                        K(n0) = K(n0) - LongInt.kmax;
                        K(n0+1) = K(n0+1) + 1;
                        n0 = n0 + 1;
                    end
                    K(i+j) = K(i+j) + d;
                    n0 = i+j;
                    while K(n0) >= LongInt.kmax
                        K(n0) = K(n0) - LongInt.kmax;
                        K(n0+1) = K(n0+1) + 1;
                        n0 = n0 + 1;
                    end
                end
            end
            while K(end) == 0
                K(end) = [];
            end
            obj3.K = K;
        end
        
        function x = divideDouble(obj1, obj2)
            
        end
    end
end

