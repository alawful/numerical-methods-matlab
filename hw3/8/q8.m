function main

% ---------- test case ----------
% System:
%   x^2 + y^2 - 1 = 0
%   x - y = 0
% Exact solution: x = y = 1/sqrt(2) ≈ 0.7071

Ffun = @(x) [x(1)^2 + x(2)^2 - 1;
             x(1) - x(2)];

Jfun = @(x) [2*x(1), 2*x(2);
             1     , -1    ];

x0   = [0.5; 0.2];
nmax = 20;
tol  = 1e-10;

[x, niter] = newtonsys(Ffun, Jfun, x0, nmax, tol);

disp('solution x =')
disp(x)
disp('iterations =')
disp(niter)

end

% ---------- Newton method for nonlinear system ----------
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
