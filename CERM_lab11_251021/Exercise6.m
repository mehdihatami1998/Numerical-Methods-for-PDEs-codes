            % Solving ODEs
            % Using ODE45 Solver

clc
clear all
close all
format long


            % Defining the problem and the Exact functions and Defaults

f = @(t, y) ( (y ^ 2) / (1 + y ^ 2) ) * exp(-10 * sin(t));

            % Defining the Default parameters
y0 = 1;
T = 4;
h = 0.01;
N = T / h;
            
            % Defining the ODE solver and options
options_ODE45 = odeset('Reltol', 10e-11,'AbsTol', 10e-12,'MaxStep', 0.01);
[tref,uref] = ode45 (f, [0 : 0.01 : T], 1, options_ODE45);


            % Solving with Heun Method
un_heun = zeros(round(N)+1,1);
tn_heun = zeros(round(N)+1,1);

un_heun(1) = y0;
tn_heun(1) = 0;

for k = 1 : N
    tn_heun(k+1,1) = tn_heun(k,1) + h;
    f1 = f(tn_heun(k,1), un_heun(k,1));
    f2 = f(tn_heun(k+1,1), un_heun(k,1)+ h * f1);
    un_heun(k+1,1) = un_heun(k,1) + 0.5 * h * (f1+f2);
end


            % Computing Absolute and Relative Errors
            % over the whole interval in the Infinity Norm

abs_err_Heun = abs(un_heun - uref);
rel_err_Heun = abs_err_Heun / abs(uref);

abs_InfNorm_Heun = norm(abs_err_Heun, inf);
rel_InfNorm_Heun = norm(rel_err_Heun, inf)




            % Solving with modified Euler (mEuler) Method

un_mEuler = zeros(round(N)+1,1);
tn_mEuler = zeros(round(N)+1,1);

un_mEuler(1) = y0;
tn_mEuler(1) = 0;

for k = 1 : N
    tn_mEuler(k+1, 1) = tn_mEuler(k, 1) + h;
    f1 = f(tn_mEuler(k, 1), un_mEuler(k, 1));
    f2 = f(tn_mEuler(k, 1) + 0.5 * h, un_mEuler(k,1) + 0.5 * h * f1);
    un_mEuler(k+1,1) = un_mEuler(k,1) + h * f2;
end


            % Computing Absolute and Relative Errors
            % over the whole interval in the Infinity Norm

abs_err_mEuler = abs(un_mEuler - uref);
rel_err_mEuler = abs_err_mEuler / abs(uref);

abs_InfNorm_mEuler = norm(abs_err_mEuler, inf);
rel_InfNorm_mEuler = norm(rel_err_mEuler, inf)




            % Solving with LeapFrog (LF) Method
un_LF = zeros(round(N)+1,1);
tn_LF = zeros(round(N)+1,1);

un_LF(1) = y0;
tn_LF(1) = 0;

un_LF(2) = uref(2);
tn_LF(2) = h;

for k = 2 : N
    tn_LF(k+1,1) = tn_LF(k, 1) + h;
    un_LF(k+1,1) = un_LF(k-1,1) + 2 * h * f(tn_LF(k), un_LF(k));
end


            % Computing Absolute and Relative Errors
            % over the whole interval in the Infinity Norm

abs_err_LF = abs(un_LF - uref);
rel_err_LF = abs_err_LF / abs(uref);

abs_InfNorm_LF = norm(abs_err_LF, inf);
rel_InfNorm_LF = norm(rel_err_LF, inf)


