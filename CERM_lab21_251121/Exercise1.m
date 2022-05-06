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
M = 1000; % whatever makes the stability condition satisfied.
dx = L / N;
dt = T / M;


% step 2
ni = @(x) 0.1 + 0.01 * cos(x);
s = @(x, t) -sin(t) * (1 + x^3/10) + (3*x^2/1000) * sin(x) * cos(x) - ...
    (0.1 + 0.01 * cos(x)) * cos(t) * 0.6 * x;

yex = @(x, t) cos(t) * (1 + 0.1 * x.^3);


c0 =@(x) 1 + 0.1 * x .^3;

% dirichlet:
g0 = @(t) cos(t);
gL = @(t) cos(t) * (1 + 0.1 * pi^3);

% Neumann BC:

% Flux BC:


% step 3
xin = [dx : dx : L-dx]'; % internal nodes for dirichlet BC
xhalf = [dx * 0.5 : dx : L - dx *0.5]'; % for ni

% step 4
t = 0
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

error1 = abs(cnew - yex(xin, T));
norm2error1 = norm(error1, 2) * dx;



%%%%%%%%% second solution



% step 1
L = pi;
T = 1;
N = 200;
M = 200000; % whatever makes the stability condition satisfied.
dx = L / N;
dt = T / M;


% step 2
ni = @(x) 0.1 + 0.01 * cos(x);
s = @(x, t) -sin(t) * (1 + x^3/10) + (3*x^2/1000) * sin(x) * cos(x) - ...
    (0.1 + 0.01 * cos(x)) * cos(t) * 0.6 * x;

yex = @(x, t) cos(t) * (1 + 0.1 * x.^3);


c0 =@(x) 1 + 0.1 * x .^3;

% dirichlet:
g0 = @(t) cos(t);
gL = @(t) cos(t) * (1 + 0.1 * pi^3);

% Neumann BC:

% Flux BC:


% step 3
xin = [dx : dx : L-dx]'; % internal nodes for dirichlet BC
xhalf = [dx * 0.5 : dx : L - dx *0.5]'; % for ni

% step 4
t = 0
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

error2 = abs(cnew - yex(xin, T));
norm2error2 = norm(error2, 2) * dx;


p_emp = -log2( norm2error2 / norm2error1)