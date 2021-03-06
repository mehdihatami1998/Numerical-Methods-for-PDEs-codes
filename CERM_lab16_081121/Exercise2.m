% CERM_Lab_081121
% Exercise 2







clear all
close all
format long

            % To define A and check if it is Diagonally Dominant
A = diag(5 * ones(200, 1), 0) + diag(-1 * ones(199, 1), 1)...
    + diag(-1 * ones(199, 1), -1) ...
    + diag(-1 * ones(198, 1), 2) ...
    + diag(-1 * ones(198, 1), -2);

A(1:5, 1:5)

m = [];

for i = 1:200
   t= sum(abs(A(i, :))) - abs(A(i, i));

   mnew = A(i, i) > t;
   m = [m; mnew];
end


if m == ones(size(m, 1))
    disp('Diagonally dominant')
else
    disp('Not Diagonally Dominant')
end


            % Compute b = Ax

Xex = ones(200, 1);

b = A * Xex;

[L, U, P] = lu(A);

y = L \ P * b;

Xlu = U \ y;


            % Compute Absolute and Relative Errors
abs_err_norm2 = norm(abs(Xex - Xlu), 2)
rel_err_norm2 = abs_err_norm2/norm(Xex, 2)

abs_err_normInf = norm(abs(Xex - Xlu), inf)
rel_err_normInf = abs_err_normInf / norm(Xex, inf)
