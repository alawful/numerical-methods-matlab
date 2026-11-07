clear
clc
close all

n = 25;
a_list = [2.0 2.2 2.4 2.6];
m = 10;

eJ = zeros(size(a_list));
eG = zeros(size(a_list));

for k = 1:length(a_list)
    a = a_list(k);

    A = a*eye(n);
    for i = 1:n-1
        A(i,i+1) = -1;
        A(i+1,i) = -1;
    end

    b = ones(n,1);
    x_true = A\b;
    x0 = zeros(n,1);

    D = diag(diag(A));
    L = tril(A,-1);
    U = triu(A,1);

    x = x0;
    for j = 1:m
        x = D\(b - (L+U)*x);
    end
    eJ(k) = norm(x - x_true, inf);

    x = x0;
    for j = 1:m
        x = (D+L)\(b - U*x);
    end
    eG(k) = norm(x - x_true, inf);
end

figure
plot(a_list, eJ, 'bo-', a_list, eG, 'r*-')
xlabel('alpha')
ylabel('final error (10 iters)')
legend('Jacobi','Gauss-Seidel')
grid on
