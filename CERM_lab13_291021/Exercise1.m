% CERM_Lab_291021
% Exercise 1
            
           
            % Solving ODEs
clc
clear all
close all
format long

f = @(t, y) -0.5 * y.^2 + 5 * y * sin(t);
y0 = 7;

T = 3;
N = 400;
h = T/N;

            % Reference Solution using ODE45
options = odeset('MaxStep',T/4000);
[tref, uref] = ode45(f,[0 : h :T], y0, options);


            % solution using Heun Method:
unHeun = zeros(N+1, 1);
tnHeun = zeros(N+1, 1);

unHeun(1, 1) = y0;
tnHeun(1, 1) = 0;

for k = 1 : N
    tnHeun(k+1, 1) = tnHeun(k, 1) + h;
    f1 = f(tnHeun(k, 1), unHeun(k, 1));
    f2 = f(tnHeun(k+1, 1), unHeun(k, 1) + h *f1);
    unHeun(k+1, 1) = unHeun(k, 1) + 0.5 * h * (f1 +f2);
end

figure(1)
plot(tref, uref, 'r', tnHeun, unHeun, 'k*')
title('Approximate solution using Heun Method')

relErrHeun = abs(uref - unHeun) / abs(uref);
normInfHeun = norm(relErrHeun, inf)

            % solution using two step Adams Bashforth method(AB2)
unAB2 = zeros( N+1, 1);
tnAB2 = zeros( N+1, 1);

unAB2(1, 1) = y0;
tnAB2(1, 1) = 0;
tnAB2(2, 1) = h;
unAB2(2, 1) = uref(2, 1);

for k = 2:N
    tnAB2(k+1, 1) = tnAB2(k, 1) + h;
    unAB2(k+1, 1) = unAB2(k, 1) + h * (1.5*f(tnAB2(k, 1),unAB2(k, 1)) ...
        - 0.5 * f(tnAB2(k-1, 1), unAB2(k-1, 1)));
end




figure(2)
plot(tref, uref, 'r', tnAB2, unAB2, 'k*')
title('Approximate solution using two step Adam Bashforth method')


relErrAB2 = abs(uref - unAB2) / abs(uref);
normInfAB2 = norm(relErrAB2, inf)


    % In case of this question, Heun method is more accurate comparing with 
    % two step Adam Bashforth method.

% why?
