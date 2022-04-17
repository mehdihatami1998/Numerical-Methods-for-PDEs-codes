% CERM_Lab_291021
% Exercise 3
% Heun Method
            % Systems of equations
            % Solve the system with Heun Method
            % Compute P_emp for Heun Method


%%%%%%%%%%%%%%%%%% Comments %%%%%%%%%%%%%%%%%%

    % Heun Method's results are totally coherent with the theory.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
format long

A = [-1 20;-30 -0.001];
g = @(t) [t .* sin(t) ; 0];
y0 = [10000000; 20000000];

f = @(t, y) A*y + g(t);
T = 2;
N = [500;1000]
            % Solving with N1 


N1 = N(1,1);
h1 = T / N1;

options = odeset('RelTol',10e-8);
[tref, uref] = ode23tb(f, [0 : h1 : T], y0, options);


unHeunN1 = zeros(N1 + 1, 2);
tnHeunN1 = zeros(N1 + 1, 1);

unHeunN1 (1,:) = y0';
tnHeunN1 (1,1) = 0;

for k = 1 : N1
    tnHeunN1(k+1, 1) = tnHeunN1(k, 1) + h1;
    f1N1 = f(tnHeunN1(k, 1)', unHeunN1(k,:)');
    f2N1 = f(tnHeunN1(k+1, 1)', unHeunN1(k,:)' + h1 *f1N1);
    unHeunN1(k+1, :) = unHeunN1(k, :)' + 0.5 * h1 * (f1N1 + f2N1);
end

plot(tref, uref, 'r', tnHeunN1, unHeunN1, 'k*')

rel_err_N1 = abs(uref - unHeunN1)/abs(uref);
normInfHeunN1 = norm(rel_err_N1, inf)


            % Solving with N2 

N2 = N(2, 1);
h1 = T / N2;

options = odeset('RelTol',10e-12);
[tref, uref] = ode23tb(f, [0 : h1 : T], y0, options);


unHeunN2 = zeros(N2 + 1, 2);
tnHeunN2 = zeros(N2 + 1, 1);

unHeunN2 (1,:) = y0';
tnHeunN2 (1,1) = 0;

for k = 1 : N2
    tnHeunN2(k+1, 1) = tnHeunN2(k, 1) + h1;
    f1N2 = f(tnHeunN2(k, 1)', unHeunN2(k,:)');
    f2N2 = f(tnHeunN2(k+1, 1)', unHeunN2(k,:)' + h1 *f1N2);
    unHeunN2(k+1, :) = unHeunN2(k, :)' + 0.5 * h1 * (f1N2 + f2N2);
end

plot(tref, uref, 'r', tnHeunN2, unHeunN2, 'k*')

rel_err_N2 = abs(uref - unHeunN2)/abs(uref);
normInfHeunN2 = norm(rel_err_N2, inf)
P_emp_Heun = -log2(normInfHeunN2/normInfHeunN1)

%%


% CERM_Lab_291021
% Exercise 3
% Crank Nicolson Method

            % Systems of equations
            % Solve the system with CN Method
            % Compute P_emp for CN Method


%%%%%%%%%%%%%%%%%% Comments %%%%%%%%%%%%%%%%%%
 
    % for N1 = 500 and N2 = 1000, the result of CN is okay
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
format long

A = [-1 20;-30 -0.001];
g = @(t) [t .* sin(t) ; 0];
y0 = [100; 200];

f = @(t, y) A*y + g(t);
T = 2;
N = [5000;10000]
            % Solving with N1 


N1 = N(1,1);
h1 = T / N1;

options = odeset('RelTol',10e-12);
[tref, uref] = ode23tb(f, [0 : h1 : T], y0, options);


unCNN1 = zeros(N1 + 1, 2);
tnCNN1 = zeros(N1 + 1, 1);

unCNN1 (1,:) = y0';
tnCNN1 (1,1) = 0;

for k = 1 : N1
    tnCNN1(k+1, 1) = tnCNN1(k, 1) + h1;

     g = @(y) y - unCNN1(k, :)' - 0.5 * h1 * (f(tnCNN1(k+1,1), y) + ...
        f(tnCNN1(k, 1), unCNN1(k, :)'));
    unCNN1(k+1, :) = fsolve(g, unCNN1(k, :)');
end



plot(tref, uref, 'r', tnCNN1, unCNN1, 'k*')

rel_err_N1 = abs(uref - unCNN1)/abs(uref);
normInfCNN1 = norm(rel_err_N1, inf)



            % Solving with N2 

N2 = N(2, 1);
h2 = T/N2;

options = odeset('RelTol',10e-12);
[tref, uref] = ode23tb(f, [0 : h2 : T], y0, options);


unCNN2 = zeros(N2 + 1, 2);
tnCNN2 = zeros(N2 + 1, 1);

unCNN2 (1,:) = y0';
tnCNN2 (1,1) = 0;

for k = 1 : N2
       tnCNN2(k+1, 1) = tnCNN2(k, 1) + h2;

     g = @(y) y - unCNN2(k, :)' - 0.5 * h2 * (f(tnCNN2(k+1,1), y) + ...
        f(tnCNN2(k, 1), unCNN2(k, :)'));
    unCNN2(k+1, :) = fsolve(g, unCNN2(k, :)', options);
end


figure(2)
plot(tref, uref, 'r', tnCNN2, unCNN2, 'k*')

rel_err_N2 = abs(uref - unCNN2)/abs(uref);
normInfCNN2 = norm(rel_err_N2, inf)
P_emp_CN = -log2(normInfCNN2/normInfCNN1)
