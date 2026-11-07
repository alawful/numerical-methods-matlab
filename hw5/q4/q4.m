h = input('h = ');
tol = input('tol = ');

x = 1.2;
f = @(t) 10*exp(-1./t).*sin(2*t);

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
        fprintf('h_best=%.15g\n', hk);
        fprintf('Dk=%.15g\n', Dk);
        break
    end

    Dkm1 = Dk;
    Dk = Dkp1;
    k = k + 1;
    hk = hk1;
end
