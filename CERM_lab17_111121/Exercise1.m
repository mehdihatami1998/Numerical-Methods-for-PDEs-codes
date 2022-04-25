% CERM_Lab_111121
% Exercise 1

clc
clear all
close all
format long


            % Definig A

A = diag(6 * ones(400,1), 0) + diag(-1 * ones(399,1), -1) + diag(-1 * ones(399,1), 1) ...
    + diag(ones(396, 1), 4) + diag(ones(396, 1), -4);

A(1:5, 1:5)


            % Check if it is diagonally dominant

m = [];
for i = 1 : 400
    t = sum(abs(A(i, :))) - abs(A(i, i));
    mnew = abs(A(i, i)) > t;
    m = [m; mnew];
end
    

if m == ones(size(m, 1))
    disp('Diagonally dominant')
else
    disp('Not diagonally dominant')
end


            % to check if it is symmetric positive definite
try
    chol(A);
    disp('it is symmetric positive definite')
catch
    disp('it is not symmetric positive definite')
end

Xex = [];
for i = 1 : 200
    Xex=[Xex; 1; -1];
end


b = A * Xex;
U = chol(A);
L = U';

y = L\b;
Xchol = U \ y;




            % Compute Absolute and Relative Errors

abs_err_norm2 = norm(abs(Xex - Xchol), 2)
rel_err_norm2 = abs_err_norm2/norm(Xex, 2)

abs_err_normInf = norm(abs(Xex - Xchol), inf)
rel_err_normInf = abs_err_normInf / norm(Xex, inf)


            % apriori and a posteriori error estimates
priori = cond(A) * eps
residual = b - A * Xchol;
posteriori = cond(A) * norm(residual)/ norm(b)

