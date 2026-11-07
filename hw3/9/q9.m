function main

Ffun = @(x) [6*x(1)^3 + x(1)*x(2) - 3*x(2)^3 - 4;
             x(1)^2   - 18*x(1)*x(2)^2 + 16*x(2)^3 + 1];

Jfun = @(x) [18*x(1)^2 + x(2),          x(1) - 9*x(2)^2;
             2*x(1) - 18*x(2)^2,   -36*x(1)*x(2) + 48*x(2)^2];

x0   = [2; 2];
nmax = 5;         % exactly 5 Newton iterations
tol  = 0;         % disable early stopping

[x, niter] = newtonsys(Ffun, Jfun, x0, nmax, tol);

disp('x =')
disp(x)
disp('niter =')
disp(niter)

end

function [x, niter] = newtonsys(Ffun, Jfun, x0, nmax, tol)

x = x0(:);
niter = 0;

for k = 1:nmax
    F = Ffun(x);
    J = Jfun(x);
    dx = -J\F;
    x = x + dx;
    niter = k;
    if norm(dx, inf) < tol
        break
    end
end

end
