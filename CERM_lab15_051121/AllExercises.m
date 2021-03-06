% CERM_Lab_041121
% Exercise 3
% second order equation





%%
            % Exercise 1
            % a
A = [50 1 3; 1 6 2;1 0 1]
B = [50 1;3 20;10 4]
C = [1 2 3;2 4 5]

D = eye(3) + B*C


            % b
A == A'

As = (A + A') / 2
Aas = (A - A') / 2

As == As'

Aas == -Aas'

            
            % c
A*D == D*A

A*D - D*A

            % d
E = eye(3) + 2 * A' * A + 3 * A ^ 3

det(E)

abs(inv(A) - A^(-1)) < 10^-9 * ones(3, 3)

inv(E)
E^(-1)


%%
            % Exercise 2
format short
A1 = [D E;-E' inv(E)] 

A2 = [eye(6) zeros(6) A1;2*eye(6) A1 eye(6);A1' zeros(6) 3*eye(6)]

size(A1) 

size (A2)

            % Exercise 3

% using toeplitz function:

r = [1 2 3 4 5];
A_toeplitz = toeplitz(r)

% using diag() function:
s = 1 * ones(5, 1);
n = 2 * ones(4, 1);
m = 3 * ones(3, 1);
o = 4 * ones(2, 1);
p = 5 * ones(1, 1);

A_diag = diag(s, 0) + diag(n, 1) + diag(n, -1) + diag(m, 2) + diag(m, -2) + ...
    diag(o, 3) + diag(o, -3) + diag(p, 4) + diag(p, -4)


% using for cycles:
A_fors = zeros(5, 5);
for i = 1:5
    
    A_fors = A_fors + diag(i * ones(6-i, 1),i-1) + diag((i+1) .* ones(5-i,1), -i)% 
    i = i-1;
end
A = A_fors;

% Extract upper and lower traingular

upperA = triu(A)
lowerA = tril(A)

triu(A, 1)
tril(A, 1)

det(upperA)
det(lowerA)

eig(A)

%%
            % Exercise 4
v = [];
for i  = 1 :50
    v = [v;2;4];
end
v



A = diag(v) + diag(1 * ones(99,1), 1) + diag(-1 * ones(99,1), -1);
A(1:5, 1:5)


B = (-3*A + 2*A^2) / (eye(100) + 4 * A - A^4);
B(1:5, 1:5)

Xex = -1 * ones(100,1);
d = B * Xex;
d(1:5)

x1 = B^(-1) * d

x2 = inv(B) * d

x3 = B\ d


%%
            % Exercise 5
% a
m = hilb(7)

n = m(:, 3)
o = ones(1,7);

m(:, 3) = o
m(5,1:3) = ones(3,1)


% b
diag(m)

% c
det(m)


% d
Xex = ones(1, 7);
b = m .* Xex;

x1 = inv(m) * b

x2 = m\b

abs(x1 - x2) <10e-7





