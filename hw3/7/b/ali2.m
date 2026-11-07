clear
close all
clc

n = 25;
a = 2.0;

A = a*eye(n);
for i = 1:n-1
    A(i,i+1) = -1;
    A(i+1,i) = -1;
end

b = ones(n,1);
x_true = A\b;

x0 = zeros(n,1);
m = 100;

D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);

x = x0;
errJ = zeros(m,1);
for k = 1:m
    x = D\(b - (L+U)*x);
    errJ(k) = norm(x - x_true, inf);
end

x = x0;
errG = zeros(m,1);
for k = 1:m
    x = (D+L)\(b - U*x);
    errG(k) = norm(x - x_true, inf);
end

it = 1:m;
figure;
semilogy(it, errJ, 'bo-', it, errG, 'r*-');
xlabel('iteration');
ylabel('error (inf-norm)');
legend('Jacobi','Gauss-Seidel');
grid on;