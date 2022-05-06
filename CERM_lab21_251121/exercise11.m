% CERM_lab_251121
% Exercise 1

clc
clear all
close all
format long


% step 1
L = pi;
T = 1;
N = 100;
dx = L / N;


% step 2
ni = @(x) 0.1 + 0.01 * cos(x);
s = @(x, t) -sin(t) * (1 + x^3/10) + (3*x^2/1000) * sin(x) * cos(x) - ...
    (0.1 + 0.01 * cos(x)) * cos(t) * 0.6 * x;

yex = @(x, t) cos(t) * (1 + 0.1 * x.^3);


c0 =@(x) 1 + 0.1 * x .^3;

% dirichlet:
g0 = @(t) cos(t);
gL = @(t) cos(t) * (1 + 0.1 * pi^3);


% Flux BC:


% step 3
xin = [dx : dx : L-dx]'; % internal nodes for dirichlet BC
xhalf = [dx * 0.5 : dx : L - dx *0.5]'; % for ni


% check if the stability condition is satisfied:

M = N;
dt = T / M;
nimax = max(ni(xhalf));
lambda = dt * nimax / dx^2;

while lambda > 0.5
    M = M +10;
    dt = T / M;
    lambda = dt * nimax / dx^2;
end
M
lambda;


% step 4 with dirichlet Boundary Conditioni
t = 0;
cold = c0(xin);
cnew = cold;

for n = 1 : M

    cnew(1) = cold(1) + (dt/dx^2) * (ni(xhalf(2)) * (cold(2) - cold(1)) ...
                                        -ni(xhalf(1)) *   (cold(1) - g0(t)))...
                                        + dt * s(xin(1), t);

    for i = 2 : N-2

        cnew(i) = cold(i) + (dt/dx^2) * (ni(xhalf(i+1)) * (cold(i + 1) - cold(i)) ...
                                        -ni(xhalf(i)) *   (cold(i) - cold(i-1)))...
                                        + dt * s(xin(i), t);
    end

      cnew(N-1) = cold(N-1) + (dt/dx^2) * (ni(xhalf(N)) * (gL(t) - cold(N-1)) ...
                                        -ni(xhalf(N-1)) *   (cold(N-1) - cold(N-1)))...
                                        + dt * s(xin(N-1), t);

    t = t+dt;
    cold = cnew;
end

errorDirichlet = abs(cnew - yex(xin, T));
norm2errorDirichlet = norm(errorDirichlet, 2) * dx


%%
% Neuman Boundary condition:


% CERM_lab_251121
% Exercise 1


clear all
close all
format long


% step 1
L = pi;
T = 1;
N = 100;
dx = L / N;


% step 2
ni = @(x) 0.1 + 0.01 * cos(x);
s = @(x, t) -sin(t) * (1 + x^3/10) + (3*x^2/1000) * sin(x) * cos(x) - ...
    (0.1 + 0.01 * cos(x)) * cos(t) * 0.6 * x;

yex = @(x, t) cos(t) * (1 + 0.1 * x.^3);


c0 =@(x) 1 + 0.1 * x .^3;

% Neuman
g0 = @(t) 0
gL = @(t) cos(t) * ( 0.1 * 3 * pi^2);





% step 3
xin = [dx : dx : L-dx]'; % internal nodes for dirichlet BC
xhalf = [dx * 0.5 : dx : L - dx *0.5]'; % for ni


% check if the stability condition is satisfied:

M = N;
dt = T / M;
nimax = max(ni(xhalf));
lambda = dt * nimax / dx^2;

while lambda > 0.5
    M = M +10;
    dt = T / M;
    lambda = dt * nimax / dx^2;
end
M
lambda;


% step 4 with Neuman Boundary Conditioni
t = 0;
cold = c0(xin);
cnew = cold;

for n = 1 : M

    cnew(1) = cold(1) + (dt/dx^2) * (ni(xhalf(2)) * (cold(2) - cold(1)) ...
                                        -ni(xhalf(1)) * (g0(t) * dx))...
                                        + dt * s(xin(1), t);

    for i = 2 : N-2

        cnew(i) = cold(i) + (dt/dx^2) * (ni(xhalf(i+1)) * (cold(i + 1) - cold(i)) ...
                                        -ni(xhalf(i)) *   (cold(i) - cold(i-1)))...
                                        + dt * s(xin(i), t);
    end

      cnew(N-1) = cold(N-1) + (dt/dx^2) * (ni(xhalf(N)) * (gL(t)*dx) ...
                                        -ni(xhalf(N-1)) *   (cold(N-1) - cold(N-1)))...
                                        + dt * s(xin(N-1), t);

    t = t+dt;
    cold = cnew;
end

errorNeuman = abs(cnew - yex(xin, T));
norm2errorNeuman = norm(errorNeuman, 2) * dx

%%
% Flux Boundary Condition

% CERM_lab_251121
% Exercise 1


clear all
close all
format long


% step 1
L = pi;
T = 1;
N = 100;
dx = L / N;


% step 2
ni = @(x) 0.1 + 0.01 * cos(x);
s = @(x, t) -sin(t) * (1 + x^3/10) + (3*x^2/1000) * sin(x) * cos(x) - ...
    (0.1 + 0.01 * cos(x)) * cos(t) * 0.6 * x;

yex = @(x, t) cos(t) * (1 + 0.1 * x.^3);


c0 =@(x) 1 + 0.1 * x .^3;



% Flux BC:

g0 = @(t) 0;
gL = @(t) cos(t) * (27 * pi ^ 2 / 1000);

% step 3
xin = [dx : dx : L-dx]'; % internal nodes for dirichlet BC
xhalf = [dx * 0.5 : dx : L - dx *0.5]'; % for ni


% check if the stability condition is satisfied:

M = N;
dt = T / M;
nimax = max(ni(xhalf));
lambda = dt * nimax / dx^2;

while lambda > 0.5
    M = M +10;
    dt = T / M;
    lambda = dt * nimax / dx^2;
end
M
lambda;


% step 4 with Flux Boundary Conditioni
t = 0;
cold = c0(xin);
cnew = cold;

for n = 1 : M

    cnew(1) = cold(1) + (dt/dx^2) * (ni(xhalf(2)) * (cold(2) - cold(1)) ...
                                        -ni(xhalf(1))) *(g0(t) * dx)...
                                        + dt * s(xin(1), t);

    for i = 2 : N-2

        cnew(i) = cold(i) + (dt/dx^2) * (ni(xhalf(i+1)) * (cold(i + 1) - cold(i)) ...
                                        -ni(xhalf(i)) *   (cold(i) - cold(i-1)))...
                                        + dt * s(xin(i), t);
    end

      cnew(N-1) = cold(N-1) + (dt/dx^2) * (ni(xhalf(N))) * (gL(t) * dx) ...
                                        -ni(xhalf(N-1)) *   (cold(N-1) - cold(N-1))...
                                        + dt * s(xin(N-1), t);

    t = t+dt;
    cold = cnew;
end

errorFlux = abs(cnew - yex(xin, T));
norm2errorFlux = norm(errorFlux, 2) * dx
