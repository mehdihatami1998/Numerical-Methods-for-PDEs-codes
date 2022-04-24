% CERM_Lab_081121
% Exercise 2

clc
clear all
close all
format long

            % check if the matrix is diagonally dominant by rows
alpha = 1;
n = 20;

A = alpha .* eye(n) + hilb(n);

m = [];
for i = 1 : n
    t = sum(abs(A(i, :))) - abs(A(i, i));
    mnew = abs(A(i, i)) > t;
    m = [m; mnew];
end
    

if m == ones(size(m, 1))
    disp('Diagonally dominant')
else
    disp('Not diagonally dominant')
end


            % Check if it is symmetric

p = [];
for i = 1 : n
    for j = 1:n
    pnew = A(i, j) == A(j, i);
    p = [p; pnew];
    end
end

if p == ones(size(p))
    disp('Symmetric')
else
    disp('asymmetric')
end


            % Solve system with LU factorization


Xex = ones(20,1);
Xex(1, 1) = -1;
Xex(20, 1) = -1;

b = A * Xex;

[L, U, P] = lu(A);
y = L \ P * b;

Xlu = U \ y;




            % Compute Absolute and Relative Errors
abs_err_norm2 = norm(abs(Xex - Xlu), 2)
rel_err_norm2 = abs_err_norm2/norm(Xex, 2)

abs_err_normInf = norm(abs(Xex - Xlu), inf)
rel_err_normInf = abs_err_normInf / norm(Xex, inf)

