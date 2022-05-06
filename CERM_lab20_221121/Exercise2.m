% CERM_lab_21121
% Exercise 2

%%%%%%%% comment:
% here we are having time disceretization of order 1 and space
% discertization of order 2, so the empirical convergence should be a
% number between these two 


clear all
close all
clc

l = @(v, dx) sqrt(sum(dx * abs(v) .^ 2));

L = 4;
T = 5;
v = 0.1;
N = 50;

dx = L / N;
xn = [0 : dx : L-dx]';


c0 = @(x) 0;
c0hat = fft(zeros(size(xn)));
fc0 = ifft(c0hat);

s = @(x) exp(-((x-L/3)/(L/20)) .^2);
shat = fft(s(xn));
fs = ifft(shat);


            % exact solution with separation of variables


kk = [0 : N/2, -N/2+1 : -1]';
omega = 2 * pi / L;

omegak = omega * kk;

alphak = v * (omegak * j).^2 ;


chat = zeros(N,1);
chat(1) =  shat(1) * T;

chat (2:N) = shat(2:N) .* (exp(alphak(2:N) * T)-1)./ alphak(2:N);


cheat = ifft(chat);

plot(xn, fs, 'b', xn, cheat, 'r')

legend('source function', 'exact answer')

hold on

            % Separation of variables

N = 50;
M = 500;
dt = T / M;
dx = L / N;



cold = zeros(N,1);
cnew = cold;
t = 0;

for n = 1 : M

    cnew(1) = cold(1) + dt * v * (cold(2) - 2*cold(1) + cold(N))./(dx.^2) + dt * s(xn(1)); %s(xper(1), t)


    for i = 2 : N-1

        cnew(i) = cold(i) + dt * v * (cold(i+1) - 2*cold(i) + cold(i-1))./ (dx.^2) + dt * s(xn(i));%s(xper(i), t)

    end

    cnew(N) = cold(N) + dt * v * (cold(1) - 2*cold(N) + cold(N-1))./ (dx.^2) + dt * s(xn(N));%s(xper(), t)

    t = t + dt;
    cold = cnew;

end
error1 = cheat - cnew;
abserror1 = abs(cheat - cnew);
norm2_1 = norm(abserror1)*sqrt(dx);
% 
% relError1Norm2 = l(cheat - cnew, dx) / l(cheat, dx);
% relError1NormInf = norm(cheat - cnew, inf)/ norm(cheat, inf);
% 


plot(xn, cnew, 'b*')
legend('source function', 'exact answer', 'Explicit Euler')




%%%%%%%%%%%%%%%%%%%% do it again with  N = 100 and M = 1000 to find p_emp

clc

l = @(v, dx) sqrt(sum(dx * abs(v) .^ 2));

L = 4;
T = 5;
v = 0.1;
N = 100;

dx = L / N;
xn = [0 : dx : L-dx]';


c0 = @(x) 0;
c0hat = fft(zeros(size(xn)));
fc0 = ifft(c0hat);

s = @(x) exp(-((x-L/3)/(L/20)) .^2);
shat = fft(s(xn));
fs = ifft(shat);


            % exact solution with separation of variables


kk = [0 : N/2, -N/2+1 : -1]';
omega = 2 * pi / L;

omegak = omega * kk;

alphak = v * (omegak * j).^2 ;


chat = zeros(N,1);
chat(1) =  shat(1) * T;

chat (2:N) = shat(2:N) .* (exp(alphak(2:N) * T)-1)./ alphak(2:N);


cheat = ifft(chat);
figure()
plot(xn, fs, 'b', xn, cheat, 'r')

legend('source function', 'exact answer')

hold on

            % Separation of variables

N = 100;
M = 1000;
dt = T / M;
dx = L / N;



cold = zeros(N,1);
cnew = cold;
t = 0;

for n = 1 : M

    cnew(1) = cold(1) + dt * v * (cold(2) - 2*cold(1) + cold(N))./(dx.^2) + dt * s(xn(1));


    for i = 2 : N-1

        cnew(i) = cold(i) + dt * v * (cold(i+1) - 2*cold(i) + cold(i-1))./ dx.^2 + dt * s(xn(i));

    end

    cnew(N) = cold(N) + dt * v * (cold(1) - 2*cold(N) + cold(N-1))./ dx.^2 + dt * s(xn(N));

    t = t + dt;
    cold = cnew;

end

error2 = cheat - cnew;
abserror2 = abs(cheat - cnew);
norm2_2 = norm(abserror2)*sqrt(dx);
% relError2Norm2 = l(cheat - cnew, dx) / l(cheat, dx);
% relError2NormInf = norm(cheat - cnew, inf)/ norm(cheat, inf);

% p_emp = -log2(l(error2, L/1000)/l(error1, L/500))
p_emp = -log2(norm2_2/norm2_1)

plot(xn, cnew, 'b*')
legend('source function', 'exact answer', 'Explicit Euler')

%pemp not working correctly