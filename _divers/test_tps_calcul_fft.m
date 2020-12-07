for k = 0:8
    n = 10^k;
    
    a = rand(1, n);
    tic;
    fft(a);
    t = toc;
    
    disp('~~~~~~~~~~~~~~');
    fprintf('n = %g\n', n);
    fprintf('T = %.2es\n', t);
    fprintf('T/n.log(n) = %E\n', t/(n*log(n)));
end


disp(' ');
disp(' ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
disp(' ');


for n = logspace(log10(1), log10(100), 7)
    n = round(n);
    
    a = rand(1e8, n);
    
    tic;
    for k = 1:size(a, 2)
        fft(a(:, k));
    end
    t = toc;
    
    tic;
    fft(a);
    T = toc;
    
    disp('~~~~~~~~~~~~~~');
    fprintf('n = %g\n', n);
    fprintf('T/nt = %.f%%\n', 100*T/t);
end