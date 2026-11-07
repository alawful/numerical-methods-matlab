fs = input('f(t)= ','s');
f = str2func(fs);

h = input('h = ');
n = input('n = ');

x = 1;

D = zeros(n+1,n+1);

for j = 0:n
    hj = h*10^(-j);
    D(j+1,1) = (f(x+hj) - f(x-hj)) / (2*hj);
end

for k = 1:n
    for j = k:n
        D(j+1,k+1) = D(j+1,k) + (D(j+1,k) - D(j,k)) / (100^k - 1);
    end
end

fprintf('D(n,n)=%.15g\n', D(n+1,n+1));
fprintf('exact=%.15g\n', cos(x));
fprintf('abs_error=%.15g\n', abs(D(n+1,n+1)-cos(x)));
