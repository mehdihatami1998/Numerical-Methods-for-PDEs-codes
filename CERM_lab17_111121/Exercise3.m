% CERM_Lab_111121
% Exercise 3

clc
clear all
close all
format short

            % Defining B

B = diag(4 * ones(300, 1), 0) + diag(ones(299, 1), 1) + diag(-1 * ones(299, 1), -1) ...
    + diag(0.5 * ones(290, 1), -10) + diag(-0.3 * ones(291, 1), 9);
B(1:11, 1:11)

C = B(1:100, 1:100);
A = C' * C;

           % to check if it is symmetric positive definite
try
    chol(A);
    disp('it is symmetric positive definite')
catch
    disp('it is not symmetric positive definite')
end


            % Define b and Xex and solving the equation
Xex = [1 : 1 : 100]';

b = A * Xex;
U = chol(A);
L = U';

y = L \ b;
Xchol = U \y;



            % Compute Absolute and Relative Errors

abs_err_norm2 = norm(abs(Xex - Xchol), 2)
rel_err_norm2 = abs_err_norm2/norm(Xex, 2)

abs_err_normInf = norm(abs(Xex - Xchol), inf)
rel_err_normInf = abs_err_normInf / norm(Xex, inf)


            % apriori and a posteriori error estimates
priori = cond(A) * eps
residual = b - A * Xchol;
posteriori = cond(A) * norm(residual)/ norm(b)


