t = 1:10;
y = [6.2 9.5 12.3 13.9 14.6 13.5 13.3 12.7 12.4 11.9];

z = log(y./t);
A = [ones(length(t),1) t'];
coef = A\z';

a = coef(1);
b = coef(2);

c1 = exp(a);
c2 = b;

disp(c1)
disp(c2)
