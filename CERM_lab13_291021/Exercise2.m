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


%%%%%%%%%%%%%%%%%% Comments %%%%%%%%%%%%%%%%%%

    % If we set N1 = 10 and N2 = 20, the P_emp would be exactly 2,
    % But if we increase N1 and N2 to about 500 and 1000, the truncation
    % error will play a role here and because the steps are very small and
    % our error is of the same order of convergance with machine accuracy,
    % So the P_emp wouldn't work in this case.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
            % Solving ODEs
clc
clear all
close all
format long

f = @(t, y) -0.5 * y.^2 + 5 * y * sin(t);
y0 = 7;


T = 3;
N1 = 10;
h1 = T/N1;

N2 = 20;
h2 = T/N2;


            % a reference solution using ode45

[tnref1, unref1] = ode45(f,[0: h1 : T], y0);

            %
            % solution using BDF3 Method with N1:
            %
unBDF3N1 = zeros(N1+1, 1);
tnBDF3N1 = zeros(N1+1, 1);

unBDF3N1(1, 1) = y0;
tnBDF3N1(1, 1) = 0;


            % Using Cranck-Nicolson method to provide first initial condition
unBDF3N1(2,1) = unref1(2,1);
tnBDF3N1(2, 1) = h1;

m = @(n) n - unBDF3N1(1, 1) - 0.5 * h1 * (f(tnBDF3N1(2, 1),n) ...
    + f(tnBDF3N1(1, 1), unBDF3N1(1, 1)));
 option = optimoptions('fsolve', 'Display', 'none', FunctionTolerance=10e-8);
unBDF3N1(2, 1) = fsolve(m, unBDF3N1(1, 1), option);


            % Using Crank-Nicolson method to provide second initial
            % condition



tnBDF3N1(3, 1) = 2 * h1;
s = @(t) t - unBDF3N1(2, 1) - 0.5 * h1 * (f(tnBDF3N1(3, 1),t) ...
    + f(tnBDF3N1(2, 1), unBDF3N1(2, 1)));

option = optimoptions('fsolve', 'Display', 'none', FunctionTolerance=10e-8);
unBDF3N1(3, 1) = fsolve(s, unBDF3N1(2, 1), option);



for k = 3 : N1
    tnBDF3N1(k+1, 1) = tnBDF3N1(k, 1) + h1;
    g = @(y) y ...
        - (18/11) * unBDF3N1(k, 1) ...
        + (9 /11) * unBDF3N1(k-1, 1) ...
        - (2 / 11)* unBDF3N1(k-2, 1) ...
        - (6 /11)  * h1 * f(tnBDF3N1(k+1, 1), y);

    unBDF3N1(k+1, 1) = fsolve(g, unBDF3N1(k, 1), option);
    
end


figure(1)
plot(tnref1, unref1, 'r', tnBDF3N1, unBDF3N1, 'k*')
title('Approximate solution using BDF3N1 Meth1od')


relErrN1 = abs(unref1 - unBDF3N1) / abs(unref1);
NormInfN1 = norm(relErrN1, inf)

            %
            % solution using BDF3 Method with N2:
            %



            % solution using BDF3N2 Method:

unBDF3N2 = zeros(N2+1, 1);
tnBDF3N2 = zeros(N2+1, 1);

unBDF3N2(1, 1) = y0;
tnBDF3N2(1, 1) = 0;


            % Using Cranck-Nicolson method to provide first initial condition
tnBDF3N2(2, 1) = h2;

m = @(n) n - unBDF3N2(1, 1) - 0.5 * h2 * (f(tnBDF3N2(2, 1),n) ...
    + f(tnBDF3N2(1, 1), unBDF3N2(1, 1)));
 option = optimoptions('fsolve', 'Display', 'none', FunctionTolerance=10e-8);
unBDF3N2(2, 1) = fsolve(m, unBDF3N2(1, 1), option);


            % Using Crank-Nicolson method to provide second initial
            % condition



tnBDF3N2(3, 1) = 2 * h2;
s = @(t) t - unBDF3N2(2, 1) - 0.5 * h2 * (f(tnBDF3N2(3, 1),t) ...
    + f(tnBDF3N2(2, 1), unBDF3N2(2, 1)));

option = optimoptions('fsolve', 'Display', 'none', FunctionTolerance=10e-8);
unBDF3N2(3, 1) = fsolve(s, unBDF3N2(2, 1), option);

[tnref2, unref2] = ode45(f,[0: h2 : T], y0);



for k = 3 : N2
    tnBDF3N2(k+1, 1) = tnBDF3N2(k, 1) + h2;
    g = @(y) y ...
        - (18/11) * unBDF3N2(k, 1) ...
        + (9 /11) * unBDF3N2(k-1, 1) ...
        - (2 / 11)* unBDF3N2(k-2, 1) ...
        - (6 /11)  * h2 * f(tnBDF3N2(k+1, 1), y);

    unBDF3N2(k+1, 1) = fsolve(g, unBDF3N2(k, 1), option);

end


figure(1)
plot(tnref2, unref2, 'r', tnBDF3N2, unBDF3N2, 'k*')
title('Approximate solution using BDF3N2 Method')

relErrN2 = abs(unref2 - unBDF3N2) / abs(unref2);
NormInfN2 = norm(relErrN2, inf)

P_emp = - log2(NormInfN2 / NormInfN1)
