N = 10;

M = zeros(N);
M(1, 1:2) = [-2 1];
M(end, end-1:end) = [1 -2];
for kl = 2:N-1
    M(kl, kl-1:kl+1) = [1 -2 1];
end

% disp(M);
disp(det(M));
% disp(inv(M));
disp(eig(M));