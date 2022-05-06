% CERM_lab_21121
% Exercise 1



%%%%%% comment %%%%%
% When computing the l2 norm of errors in this section, we should multiply 
% the norm by *sqrt(dx) 
% norm(err,2)


%%%%%%% comment #2:
% here we are using 2 different discretization, one in time, which is first
% order, and then in space, which is 2nd order, so the order of convergance
% for us should be 1, accodring to the order of disceritzation in time,
%BUT, if we make the time step very small, say 0.05, it's error would be so
%small that we can say it's negligible, and the order of convergence that
%we get for this function would be the order of convergence for the space
%discretization which is 2.


            % Separation of variables
            % Periodic boundary condition
clc
clear all
close all
format long
 
        
            % Norm 2 definition
l = @(v, dx) sqrt(sum(dx * abs(v) .^ 2));

L = 10;
T = 5;
v = 0.05;


N = 100;
dx = L / N;
xper = [0 : dx : L-dx]';

            % Initial Datum
c0 = @(x) 10 * exp(-((x-L/2)/(L/10)).^2);
c0hat = fft(c0(xper));
fc0 = ifft(c0hat);

kk = [0 : N/2, -N/2 + 1 : -1]';
omegak = (2 * pi / L) * kk;
alphak = v * (j * omegak).^2;

chat = c0hat .* exp(alphak * T);
Cinit = ifft(chat);

            % Plotting the figures

figure()
plot(xper, fc0, 'b', xper, Cinit, 'kd')

hold on
            % Numerical solution by Explicit Euler method

M = 200;
dt = T / M;

            % check stability

if dt * v / (dx.^2) < 0.5
    display ('Stable Condition')
else
    display ('Unstable Conditon')
end


t = 0;
cold = c0(xper);
cnew = cold;

for n = 1:M

    cnew(1) = cold(1) + v * dt * (cold(2) - 2 * cold(1) + cold(N)) / dx.^2;

    for i = 2 : N-1

        cnew(i) = cold(i) + v * dt * (cold(i+1) - 2*cold(i) + cold(i-1))/ dx.^2;

    end

    cnew(N) = cold(N) + v * dt * (cold(1) - 2 * cold(N) + cold(N-1)) / dx.^2;
    cold = cnew;
    t = t + dt;
end

plot(xper, cnew,'r-')
legend('Initial moment', 'Separation of Variables', 'Explicit Euler')



error = Cinit - cnew;
absErr = l(error, dx)
relErr = norm(error) / norm(Cinit)

figure()
plot(xper, error, 'ro')


            % Repeat calculations with N = 200, M = 400


N = 200;
dx = L / N;
xper = [0 : dx : L-dx]';

            % Initial Datum
c0 = @(x) 10 * exp(-((x-L/2)/(L/10)).^2);
c0hat = fft(c0(xper));
fc0 = ifft(c0hat);

kk = [0 : N/2, -N/2 + 1 : -1]';
omegak = (2 * pi / L) * kk;
alphak = v * (j * omegak).^2;

chat = c0hat .* exp(alphak * T);
Cinit = ifft(chat);

            % Plotting the figures

figure()
plot(xper, fc0, 'b', xper, Cinit, 'kd')

hold on
            % Numerical solution by Explicit Euler method

M = 400;
dt = T / M;

            % check stability

if dt * v / (dx.^2) < 0.5
    display ('Stable Condition')
else
    display ('Unstable Conditon')
end


t = 0;
cold = c0(xper);
cnew = cold;

for n = 1:M

    cnew(1) = cold(1) + v * dt * (cold(2) - 2 * cold(1) + cold(N)) / dx.^2;

    for i = 2 : N-1

        cnew(i) = cold(i) + v * dt * (cold(i+1) - 2*cold(i) + cold(i-1))/ dx.^2;

    end

    cnew(N) = cold(N) + v * dt * (cold(1) - 2 * cold(N) + cold(N-1)) / dx.^2;
    cold = cnew;
    t = t + dt;
end

plot(xper, cnew,'r-')
legend('Initial moment', 'Separation of Variables', 'Explicit Euler')



error1 = Cinit- cnew;
absErr1 = l(error1, dx)

relErr = norm(error1) / norm(Cinit)

p_emp= -log2(absErr1/absErr)
figure()
plot(xper, error1, 'ro')
