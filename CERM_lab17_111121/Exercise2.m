% CERM_Lab_081121
% Exercise 4

clc
clear all
close all
format long


%%%%%%%%%%% after n = 14 the matrix won't be positive definite  %%%%%%%%%
%%%%%%%%%%% after n = 11 the accuracy drops significantly %%%%%%%%%
prior = [];
for n = 1 : 11   
    A = hilb(n);
    i

    
                % Check if it is diagonally dominant
    
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
    
    
                % to check if it is symmetric positive definite
    try
        chol(A);
        disp('it is symmetric positive definite')
    catch
        disp('it is not symmetric positive definite')
    end


    Xex = [];
    for i = 1 :n
        Xex = [Xex;2];
    end

    b = A * Xex;

    U = chol(A);
    L = U';
    y = L \ b;

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

    prior = [prior;priori];
    

end
semilogy([1 : 1 :11], prior )