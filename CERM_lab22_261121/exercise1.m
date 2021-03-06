% CERM_lab_261121
% Exercise 1

clc
clear all
close all
format long
% separation of variables:
L = 10;
N = 100;
dx = L / N;

xper = [0 : dx : L-dx]';


T = 5;
ni = 0.1;

c0 = @(x) 10 * exp ( -((x-L/2)/(L/10)).^2);
c0hat = fft(c0(xper));
fc0 = ifft(c0hat);


kk = [0 : N/2, -N/2 + 1: -1]';
omega = 2 * pi / L;
omegak = omega * kk;

a2 = ni;
alphak = a2 * (omegak * j).^2;

chat = c0hat .* exp(alphak * T);
fc = ifft(chat);

% explicit method:

dt = 2/10 * dx.^2 / ni;
M = (T / dt);

dtdx2 = dt/dx.^2;
cold2 = c0(xper);
cnew2 = cold2;
t = 0;

for n = 1 : M

    cnew2(1) = cold2(1) + dtdx2  * ni *(cold2(2) -2*cold2(1) + cold2(N));

    for i = 2 : N-1

        cnew2(i) = cold2(i) + dtdx2 *  ni *(cold2(i+1) -2*cold2(i) + cold2(i-1));

    end

    cnew2(N) = cold2(N) + dtdx2 *  ni *(cold2(i+1) -2*cold2(i) + cold2(i-1)) ;

    t = t + dt;
    cold2 = cnew2;

end


figure(2)
plot(xper, fc0, xper, fc, 'ro', xper, cnew2, 'g--d')
legend('initial datum', 'separation of values', 'explicit')


errorexplicit = abs(cnew2 - fc);
normexplicit = norm(errorexplicit, inf)

% implicit Euler method M = 20 timesteps



t = 0;
cold = c0(xper);
cnew1 = cold;
M=20;
dt = T/M;
dtdx2 = dt/dx.^2;

for n = 1 : M
    b = cold;
    A= diag(ones(N,1) + 2 * ni * dtdx2 ,0)...
       -diag(ones(N-1,1) * ni * dtdx2  , 1)...
       -diag(ones(N-1,1) * ni * dtdx2  , -1);

    %because it's periodic boundary condition:

    A(1, N) = -dtdx2 * ni;
    A(N, 1) = -dtdx2 * ni;

    cnew1 = A\b;
    cold = cnew1;
    t = t + dt;
end

errorimplicit= abs(cnew1 - fc);
normimplicit = norm(errorimplicit, inf)

error20 = abs(cnew1 - fc);
secnorm20 = norm(error20, 2) * sqrt(dx);
infnorm20 = norm(error20, inf)


figure(1)
plot(xper, fc0, 'r', xper, fc, 'kd', xper, cnew1, 'b--o')
legend('Initial datum', 'separation of variables', 'implicit euler')



%repeat with N = 200, M = 40


L = 10;
N = 200;
dx = L / N;

xper = [0 : dx : L-dx]';


T = 5;
ni = 0.1;

c0 = @(x) 10 * exp ( -((x-L/2)/(L/10)).^2);
c0hat = fft(c0(xper));
fc0 = ifft(c0hat);


kk = [0 : N/2, -N/2 + 1: -1]';
omega = 2 * pi / L;
omegak = omega * kk;

a2 = ni;
alphak = a2 * (omegak * j).^2;

chat = c0hat .* exp(alphak * T);
fc = ifft(chat);



% implicit Euler method M = 20 timesteps
M = 40;
dt = T / M;



t = 0;
cold = c0(xper);
cnew3 = cold;
dtdx2 = dt/dx.^2;

for n = 1 : M
    b = cold;
    A= diag(ones(N,1) + 2 * ni * dtdx2 ,0)...
       -diag(ones(N-1,1) * ni * dtdx2  , 1)...
       -diag(ones(N-1,1) * ni * dtdx2  , -1);

    %because it's periodic boundary condition:

    A(1, N) = -dtdx2 * ni;
    A(N, 1) = -dtdx2 * ni;

    cnew3 = A\b;
    cold = cnew3;
    t = t + dt;
end

error40 = abs(cnew3 - fc);
secnorm40 = norm(error40, 2) * sqrt(dx);
infnorm40 = norm(error40, inf)




p_emp = -log2(secnorm40/secnorm20)


