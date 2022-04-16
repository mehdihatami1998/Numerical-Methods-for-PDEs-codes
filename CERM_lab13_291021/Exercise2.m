% CERM_Lab_291021
% Exercise 2
% Part a           
           
            % Solving ODEs
clc
clear all
close all
format long

f = @(t, y) -0.5 * y.^2 + 5 * y * sin(t);
y0 = 7;


T = 3;
N = 40;
h = T/N;


            % a reference solution using ode45

[tnref, unref] = ode45(f,[0: h : T], y0);


            % solution using BDF2 Method:

unBDF2 = zeros(N+1, 1);
tnBDF2 = zeros(N+1, 1);

unBDF2(1, 1) = y0;
tnBDF2(1, 1) = 0;


            % Using Cranck-Nicolson method to provide initial condition


            
tnBDF2(2, 1) = h;
m = @(n) n - unBDF2(1, 1) - 0.5 * h * (f(tnBDF2(2, 1),n) ...
    + f(tnBDF2(1, 1), unBDF2(1, 1)));
option = optimoptions('fsolve', 'Display', 'none', FunctionTolerance=10e-8);
unBDF2(2, 1) = fsolve(m, unBDF2(1, 1), option);



for k = 2 : N
    tnBDF2(k+1, 1) = tnBDF2(k, 1) + h;
    g = @(y) y - (4/3) * unBDF2(k, 1) + (1/3) * unBDF2(k-1, 1) ...
        - (2/3) * h * f(tnBDF2(k+1, 1), y);
    unBDF2(k+1, 1) = fsolve(g, unBDF2(k, 1), option);
    
    figure(1)
    plot(tnref, unref, 'r', tnBDF2, unBDF2, 'k*')
    title('Approximate solution using BDF2 Method')
    pause (0.1)

end





%%
% CERM_Lab_291021
% Exercise 2
% Part b

      
            % Solving ODEs
clc
clear all
close all
format long

f = @(t, y) -0.5 * y.^2 + 5 * y * sin(t);
y0 = 7;


T = 3;
N = 100;
h = T/N;


            % a reference solution using ode45

[tnref, unref] = ode45(f,[0: h : T], y0);


            % solution using BDF3 Method:

unBDF3 = zeros(N+1, 1);
tnBDF3 = zeros(N+1, 1);

unBDF3(1, 1) = y0;
tnBDF3(1, 1) = 0;


            % Using Cranck-Nicolson method to provide first initial condition
unBDF3(2,1) = unref(2,1);
tnBDF3(2, 1) = h;

m = @(n) n - unBDF3(1, 1) - 0.5 * h * (f(tnBDF3(2, 1),n) ...
    + f(tnBDF3(1, 1), unBDF3(1, 1)));
 option = optimoptions('fsolve', 'Display', 'none', FunctionTolerance=10e-8);
unBDF3(2, 1) = fsolve(m, unBDF3(1, 1), option);


            % Using Crank-Nicolson method to provide second initial
            % condition



tnBDF3(3, 1) = 2 * h;
s = @(t) t - unBDF3(2, 1) - 0.5 * h * (f(tnBDF3(3, 1),t) ...
    + f(tnBDF3(2, 1), unBDF3(2, 1)));

option = optimoptions('fsolve', 'Display', 'none', FunctionTolerance=10e-8);
unBDF3(3, 1) = fsolve(s, unBDF3(2, 1), option);



for k = 3 : N
    tnBDF3(k+1, 1) = tnBDF3(k, 1) + h;
    g = @(y) y ...
        - (18/11) * unBDF3(k, 1) ...
        + (9 /11) * unBDF3(k-1, 1) ...
        - (2 / 11)* unBDF3(k-2, 1) ...
        - (6 /11)  * h * f(tnBDF3(k+1, 1), y);

    unBDF3(k+1, 1) = fsolve(g, unBDF3(k, 1), option);

    plot(tnref, unref, 'r', tnBDF3, unBDF3, 'k*')
    pause (0.1)
end


figure(1)
plot(tnref, unref, 'r', tnBDF3, unBDF3, 'k*')
title('Approximate solution using BDF3 Method')


%%
% CERM_Lab_291021
% Exercise 2
% Part c


      
            % Solving ODEs
clc
clear all
close all
format long

f = @(t, y) -0.5 * y.^2 + 5 * y * sin(t);
y0 = 7;


T = 3;
N = 50000;
h = T/N;


            % a reference solution using ode45

[tnref, unref] = ode45(f,[0: h : T], y0);


            % solution using BDF3 Method:

unBDF3 = zeros(N+1, 1);
tnBDF3 = zeros(N+1, 1);

unBDF3(1, 1) = y0;
tnBDF3(1, 1) = 0;


            % Using Cranck-Nicolson method to provide first initial condition
unBDF3(2,1) = unref(2,1);
tnBDF3(2, 1) = h;

m = @(n) n - unBDF3(1, 1) - 0.5 * h * (f(tnBDF3(2, 1),n) ...
    + f(tnBDF3(1, 1), unBDF3(1, 1)));
 option = optimoptions('fsolve', 'Display', 'none', FunctionTolerance=10e-8);
unBDF3(2, 1) = fsolve(m, unBDF3(1, 1), option);


            % Using Crank-Nicolson method to provide second initial
            % condition



tnBDF3(3, 1) = 2 * h;
s = @(t) t - unBDF3(2, 1) - 0.5 * h * (f(tnBDF3(3, 1),t) ...
    + f(tnBDF3(2, 1), unBDF3(2, 1)));

option = optimoptions('fsolve', 'Display', 'none', FunctionTolerance=10e-8);
unBDF3(3, 1) = fsolve(s, unBDF3(2, 1), option);



for k = 3 : N
    tnBDF3(k+1, 1) = tnBDF3(k, 1) + h;
    g = @(y) y ...
        - (18/11) * unBDF3(k, 1) ...
        + (9 /11) * unBDF3(k-1, 1) ...
        - (2 / 11)* unBDF3(k-2, 1) ...
        - (6 /11)  * h * f(tnBDF3(k+1, 1), y);

    unBDF3(k+1, 1) = fsolve(g, unBDF3(k, 1), option);

    plot(tnref, unref, 'r', tnBDF3, unBDF3, 'k*')
    pause (0.00001)
end


figure(1)
plot(tnref, unref, 'r', tnBDF3, unBDF3, 'k*')
title('Approximate solution using BDF3 Method')

