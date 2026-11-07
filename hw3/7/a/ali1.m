n = 25;
alpha = 1.8;

A = diag(alpha*ones(n,1)) + diag(-ones(n-1,1),1) + diag(-ones(n-1,1),-1);
b = ones(n,1);

x_exact = A\b;

x0 = zeros(n,1);
nmax = 75;

D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);

Tj = -D\(L+U);
cj = D\b;
x = x0;
errJ = zeros(nmax,1);
for k = 1:nmax
    x = Tj*x + cj;
    errJ(k) = norm(x - x_exact, inf);
end

DL = D + L;
Tg = -DL\U;
cg = DL\b;
x = x0;
errG = zeros(nmax,1);
for k = 1:nmax
    x = Tg*x + cg;
    errG(k) = norm(x - x_exact, inf);
end

k = 1:nmax;
semilogy(k, errJ, 'b-o', k, errG, 'r-*');
xlabel('iteration');
ylabel('error (inf-norm)');
legend('Jacobi','Gauss-Seidel');


function [x, niter] = jgsmeth(A, b, x0, nmax, tol, P)

A = double(A);
b = double(b);
x = x0(:);

D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);

if P == 'J'
    T = -D\(L+U);
    c = D\b;
    for k = 1:nmax
        xnew = T*x + c;
        if norm(xnew - x, inf) < tol
            x = xnew;
            niter = k;
            return
        end
        x = xnew;
    end
    niter = nmax;

elseif P == 'G'
    DL = D + L;
    T = -DL\U;
    c = DL\b;
    for k = 1:nmax
        xnew = T*x + c;
        if norm(xnew - x, inf) < tol
            x = xnew;
            niter = k;
            return
        end
        x = xnew;
    end
    niter = nmax;

else
    x = x0;
    niter = 0;
end

end
