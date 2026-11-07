function main

A = [4 -1 -1 0;
    -1 4  0 1;
    -1 0  4 -1;
     0 -1 -1 4];

b = [1; 2; 0; 1];
x0 = zeros(4,1);
nmax = 100;
tol = 1e-8;

[xJ, itJ] = jgsmeth(A, b, x0, nmax, tol, 'J');
[xG, itG] = jgsmeth(A, b, x0, nmax, tol, 'G');

disp('Jacobi solution:')
disp(xJ)
disp('Jacobi iterations:')
disp(itJ)

disp('GaussSeidel solution:')
disp(xG)
disp('GaussSeidel iterations:')
disp(itG)

end

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
