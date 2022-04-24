% CERM_Lab_081121
% Exercise 4


clear all
close all
format long

v = [0.0001, 2, 3, 4, 5, 6, 7, 8, 9, 10]';

A = vander(v);


Xex = [10 : -1 : 1]';

b = A * Xex;

[L, U, P] = lu(A);

y = L \ P * b;
Xlu = U \ y;



            % Compute Absolute and Relative Errors
abs_err_norm2 = norm(abs(Xex - Xlu), 2)
rel_err_norm2 = abs_err_norm2/norm(Xex, 2)

abs_err_normInf = norm(abs(Xex - Xlu), inf)
rel_err_normInf = abs_err_normInf / norm(Xex, inf)


abs(det(A))

%%%%%%%%%%%%%%%%%% Comments %%%%%%%%%%%%%%%%%%

    % Here we computed X and compared that with Xex which was given by the
    % question and the error was about 10 ^ - 6 which is about 10 times
    % larger than the machine epsilone (10^-16), so this cannot be made by
    % the round off error or something. so there is something wrong here.
    % So we check the determinent and we see it's very large, then what?
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
