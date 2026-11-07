f = @sin;
x = 1;
h = 0.1;
tol = 1e-8;

maxk = 50;

k = 0;
hk = h*10^(-k);
Dkm1 = (f(x+hk) - f(x-hk)) / (2*hk);

k = 1;
hk = h*10^(-k);
Dk = (f(x+hk) - f(x-hk)) / (2*hk);

while k < maxk
    hk1 = h*10^(-(k+1));
    Dkp1 = (f(x+hk1) - f(x-hk1)) / (2*hk1);

    if abs(Dk - Dkm1) <= abs(Dkp1 - Dk) || abs(Dk - Dkm1) < tol
        fprintf('k=%d\n', k);
        fprintf('Dk=%.15g\n', Dk);
        fprintf('exact=%.15g\n', cos(x));
        fprintf('abs_error=%.15g\n', abs(Dk - cos(x)));
        break
    end

    Dkm1 = Dk;
    Dk = Dkp1;
    k = k + 1;
end
