% CERM_Lab_041121
% Exercise 1





%%%%%%%%%%%%%%%%%% Comments %%%%%%%%%%%%%%%%%%

    % explicit  Method is much faster than the Implicit methods generally, BDF3
    % in this specific case. but when the time steps are very large (very
    % small Ns) the explicit methods won't work, as it is the case here
    % with N = 20, but if we increase it, it would be working properly.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Systems of equations
            % Solve the system with BDF3 Method

clear all
close all
format long
A = [-30 0 -28;-29 -1 -29;0 0 -2];
g = [2;2;2];
f = @(t, y) A*y + g

y0 = [1;2;-10];
T = 4
N = 20
h = T / N;

tnBDF3 = zeros(N+1, 1);
unBDF3 = zeros(N+1, 3);

unBDF3(1,:) = y0;
tnBDF3(1,1) = 0;
options = odeset('RelTol',10^-9);
[tref, uref] = ode45 (f, [0 : h : T], y0, options);
plot(tref, uref)

tnBDF3(2, 1) = h;
unBDF3(2, :) = uref(2, :);

tnBDF3(3, 1) = 2 * h;
unBDF3(3, :) = uref(3, :);
tic
for k = 3 : N
    tnBDF3(k+1, 1) = tnBDF3(k, 1) + h;
    g = @(y) y ...
        - (18 / 11) * unBDF3(k, :)'...
        + (9 / 11)  * unBDF3(k-1, :)'...
        - (2 / 11)  * unBDF3(k-2, :)'...
        - (6 / 11)  * h * f(tnBDF3(k+1, 1), y);
    unBDF3(k+1, :) = fsolve(g, unBDF3(k, :)');


end
toc
plot(tref, uref, tnBDF3, unBDF3, 'k*')
title('BDF3 on the reference answer')


relErr1 = abs(uref(:, 1) - unBDF3(:, 1)) / abs(uref(:, 1));
normInf1 = norm(relErr1, inf);


relErr2 = abs(uref(:, 2) - unBDF3(:, 2)) / abs(uref(:, 2));
normInf2 = norm(relErr2, inf);


relErr3 = abs(uref(:, 3) - unBDF3(:, 3)) / abs(uref(:, 3));
normInf3 = norm(relErr3, inf);

NormInf = [normInf1 normInf2 normInf3]

%%
            % Use three stage Runge Kutta method (RK3)


% CERM_Lab_041121
% Exercise 3
% Heun Method
            % Systems of equations
            % Solve the system with RK3 Method

clear all
close all
format long
A = [-30 0 -28;-29 -1 -29;0 0 -2];
g = [2;2;2];
f = @(t, y) A*y + g

y0 = [1;2;-10];
T = 4
N = 20
h = T / N;

tnRK3 = zeros(N+1, 1);
unRK3 = zeros(N+1, 3);

unRK3(1,:) = y0;
tnRK3(1,1) = 0;
options = odeset('RelTol',10^-9);
[tref, uref] = ode45 (f, [0 : h : T], y0, options);

tic
for k = 1 : N
    tnRK3(k+1, 1) = tnRK3(k, 1) + h;
    f1 = f(tnRK3(k, 1), unRK3(k, :)');
    f2 = f(tnRK3(k, 1) + h/2, unRK3(k,:)' + h/2 *f1);
    f3 = f(tnRK3(k, 1) + h/2, unRK3(k,:)' + h/2 *f2);
    f4 = f(tnRK3(k+1, 1) , unRK3(k,:)' + h * f3);
    unRK3(k+1, :) = unRK3(k, :)' + h/6 *(f1 + 2*f2 +2*f3 +f4) ;
end
toc

figure(2)
plot(tref, uref, tnRK3, unRK3, 'k*')
title('Runge Kutta on the reference answers')


relErr1 = abs(uref(:, 1) - unRK3(:, 1)) / abs(uref(:, 1));
normInf1 = norm(relErr1, inf);


relErr2 = abs(uref(:, 2) - unRK3(:, 2)) / abs(uref(:, 2));
normInf2 = norm(relErr2, inf);


relErr3 = abs(uref(:, 3) - unRK3(:, 3)) / abs(uref(:, 3));
normInf3 = norm(relErr3, inf);

NormInf = [normInf1 normInf2 normInf3]
