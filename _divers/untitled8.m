N = 65;
M = N - 36;
n = 35;

% N = 65 - 20;
% M = N - 36;
% n = 35 - 10;

Pk = @(k) nchoosek(n-1, k) * nchoosek(N-n, M-k) / nchoosek(N-1, M);
Pn2 = @(n2) Pk(n-n2);

P = zeros(1, n);
for indn = 1:n
    try
        P(indn) = Pn2(indn);
    catch
        P(indn) = nan;
    end
end

stem(P);