clear; clc; close all;

f = @(x,y) x.^2 - y;

x0 = 0; xf = 4; y0 = 1;

hs = [0.25, 0.125, 0.0625];
x_table = x0:0.25:xf;

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t_ref, y_ref] = ode45(@(x,y) f(x,y), x_table, y0, opts);
y_ref = y_ref(:);

x_dense = linspace(x0, xf, 600);
[t_dense, y_dense] = ode45(@(x,y) f(x,y), x_dense, y0, opts);

figure; hold on; grid on;
plot(t_dense, y_dense, 'k-', 'LineWidth', 2);

for m = 1:length(hs)
    h = hs(m);
    N = round((xf - x0)/h);

    x = x0 + (0:N)*h;
    y = zeros(size(x));
    y(1) = y0;

    for n = 1:N
        y(n+1) = y(n) + h * f(x(n), y(n));
    end

    idx = round((x_table - x0)/h) + 1;
    y_e_tab = y(idx).';
    err_tab = abs(y_e_tab - y_ref);

    fprintf('\n=============================================\n');
    fprintf('Euler method results for h = %.6f\n', h);
    fprintf('   x        y_Euler           y_ode45           abs_error\n');
    fprintf('-------------------------------------------------------------\n');
    for k = 1:length(x_table)
        fprintf('%5.2f   %14.10f   %14.10f   %14.10e\n', ...
            x_table(k), y_e_tab(k), y_ref(k), err_tab(k));
    end

    if m == 1
        plot(x, y, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
    elseif m == 2
        plot(x, y, 's-', 'LineWidth', 1.2, 'MarkerSize', 4);
    else
        plot(x, y, 'd-', 'LineWidth', 1.2, 'MarkerSize', 4);
    end
end

xlabel('x');
ylabel('y(x)');
title('Euler solutions vs ode45 reference');
legend('ode45 (reference)', 'Euler h=0.25', 'Euler h=0.125', 'Euler h=0.0625', 'Location','northeast');
