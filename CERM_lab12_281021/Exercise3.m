% CERM_Lab_281021
% Exercise 2
            
           
            % Solving ODEs
            % Reference Solution using ODE15s
clc
clear all
close all
format long

f = @(t, y) -50 * y + 4 * y^2 + t;
y0 = 2;

T = 1;

optionsODE = odeset('RelTol',10e-10, 'MaxStep', 10e-4);
[tref, uref] = ode15s(f,[0 : 0.01 : 1], y0, optionsODE);

            % Explit Euler with N=15
tic
N = 100;
h = T / N;
un = zeros(N+1, 1);
tn = zeros(N+1, 1);

un(1,1) = y0;
tn(1,1) = 0;

for k = 1 : N
    tn(k+1, 1) = tn(k, 1) + h;
    un(k+1, 1) = un(k, 1) + h * f(tn(k, 1), un(k, 1));
end

figure(1)
plot(tref, uref,'r', tn, un, 'k*')
title('Explicit Euler method on the Reference answer')
toc
%%

    % We see that with N=15 the explicit Euler method doesn't work properly
    % but if we increase N to about 100, it will work according the theroy


                % Implicit Euler Method
tic
Ni = 100;
hi = T/Ni;

tni = zeros(Ni+1, 1);
uni = zeros(Ni+1, 1);

tni(1, 1) = 0;
uni(1, 1) = y0;

optionsImplicit = optimoptions('fsolve', 'Display','none', FunctionTolerance=10e-6)

for k = 1 : Ni
    tni(k+1, 1) = tni(k, 1) + hi;
    g = @(y) y - uni(k, 1) - h * f(tni(k+1), y);
    uni(k+1, 1) = fsolve (g, uni(k, 1), optionsImplicit);
end

toc

optionsODE1 = odeset('RelTol',10e-10, 'MaxStep', 10e-4);
[tref2, uref2] = ode15s(f,[0 : 0.01 : 1], y0, optionsODE1);


figure (2)
plot(tref2, uref2, 'r', tni, uni, 'k*')
title('Euler implicit on reference answer')

absErr = abs(tref2 - tni);
Norm2abs = norm(absErr,inf)

    % In terms of the speed, because the in the Implicit Method MATLAB
    % needs to solve an extra equation in each step, it's almost 4 times
    % slover than the explicit method.
    % But unlike the Explicit Method, for implicit methods we dont' need to
    % satisfy the equation stability and with any small number of steps the
    % method works properly.