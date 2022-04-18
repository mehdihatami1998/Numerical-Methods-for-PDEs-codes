% CERM_Lab_041121
% Exercise 3
% second order equation




clc
clear all
close all
format long 


mu = 2;
f = @(t, y) [y(2);-y(1)+ mu * (1 - y(1)^2) * y(2)];

y0 = [1;-1];
T = 20;
N = 500;
h = T/N;

options = odeset('RelTol', 10^-10, 'MaxStep', h);
tnref = [0 : h : T]';
[tnref, unref] = ode15s(f, tnref, y0, options);

figure(1)
plot(tnref, unref, 'k')
title('reference answer using ode15s')
hold on

            % Three stage Runge-Kutta method

tnRK3 = zeros(N+1, 1);
unRK3 = zeros(N+1, 2);

unRK3(1, :) = y0';
tnRK3(1, 1) = 0;

for k = 1 : N
    tnRK3(k+1, 1) = tnRK3(k, 1) + h;
    f1 = f(tnRK3(k, 1), unRK3(k,:)');

    f2 = f(tnRK3(k, 1) + h, unRK3(k, :)' + h *f1);

    f3 = f(tnRK3(k, 1) + 0.5 * h, unRK3(k, :)' + 0.25 * h * (f1+f2));

    unRK3(k+1, :) = [unRK3(k, :)' + h/6 * (f1 + f2 + 4*f3)];
end

plot(tnRK3, unRK3, 'r--')
title(' reference answers in black line and 3 stage Runge kuttah in red dashes ')

relErrRK3 = abs(unref - unRK3)/ abs(unref);
normInfRK3 = norm(relErrRK3, inf)

            % Three-step Adamas Bashforth method

tnAB3 = zeros(N+1, 1);
unAB3 = zeros(N+1, 2);

tnAB3(1, 1) = 0;
unAB3(1, :) = y0';

tnAB3(2, 1) = h;
unAB3(2, :) = unRK3(2, :);

tnAB3(3, 1) = 2 * h;
unAB3(3, :) = unRK3(3, :);

for k = 3 : N
    tnAB3 (k+1, 1) = tnAB3(k, 1) + h;
    
    unAB3(k+1, :) = [unAB3(k, :)' + h * ( ...
        (23/12) * f(tnAB3(k, 1), unAB3(k, :)') ...
        -(16/12) * f(tnAB3(k-1, 1), unAB3(k-1, :)') ...
        +(5/12) * f(tnAB3(k-2, 1), unAB3(k-2, :)'))]';
end

figure(2)
plot(tnref, unref, 'k', tnAB3, unAB3, 'r--')
title('Refernence answer using ode15s in black lines and answer using 3 step Adams Bashforth in red dashes')


relErrAB3 = abs(unref - unAB3)/ abs(unref);
normInfAB3 = norm(relErrAB3, inf)

%%
            % Section 2 : Using Cranck Nicolson and BDF2 Methods 
            % with N = 100 intervals and mu = 20

clc
clear all
close all
format long 


mu = 5;
f = @(t, y) [y(2);-y(1)+ mu * (1 - y(1)^2) * y(2)];

y0 = [1;-1];
T = 20;
N = 100;
h = T/N;

options = odeset('RelTol', 10^-10, 'MaxStep', h);
tnref = [0 : h : T]';
[tnref, unref] = ode15s(f, tnref, y0, options);


            % Crank Nicolson Method

tnCN = zeros(N+1, 1);
unCN = zeros(N+1, 2);

tnCN(1, 1) = 0;
unCN(1, :) = y0';
options = optimoptions ('fsolve', 'Display','none','FunctionTolerance',10^-7);

for k = 1 : N
    tnCN(k+1, 1) = tnref(k, 1) + h;

    g = @(y) y - unCN(k, :)' - 0.5 * h * (f(tnCN(k+1, 1), y) + ...
        f(tnCN(k, 1), unCN(k, :)'));
    unCN(k+1, :) = [fsolve(g, unCN(k, :)', options)]';
end

figure (3)
plot(tnref, unref, 'k', tnCN, unCN, 'r--')
title(' reference solution using ODE15s in red and Crank Nicolson method in red dashes')
            % for mu = 2 it works properly but for mu = 20 the function is
            % too stiff to be solved with this solvers


                        % BDF2 Methood

tnBDF2 = zeros(N+1, 1);
unBDF2 = zeros(N+1, 2);

tnBDF2(1, 1) = 0;
unBDF2(1, :) = y0';

tnBDF2(2, 1) = h;
unBDF2(2, :) = unCN(2, :);

for k = 2 : N
    tnBDF2(k+1, 1) = tnBDF2(k, 1) + h;

    g = @(y) y - (4/3) * unBDF2(k, :)' + (1/3) * unBDF2(k-1, :)' - ...
          (2/3) * h * f(tnBDF2(k+1, 1), y);

    unBDF2(k+1, :) = [fsolve(g, unBDF2(k, :)', options)]';
end


figure (4)
plot(tnref, unref, 'k', tnBDF2, unBDF2, 'r--')
title(' reference solution using ODE15s in red and BDF2 method in red dashes')

