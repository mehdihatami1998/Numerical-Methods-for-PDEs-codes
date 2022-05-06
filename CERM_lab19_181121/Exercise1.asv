% CERM_lab_181121
% Exercise 1

clc
clear all
close all
format long
            
            
            % Initial condition
L = 10;
T = 5;
N = 200;
nmodes = 100;
h = L / N;
xn = [0 : h : L-h]';

c01 = @(x) exp(-((x-L/2) / (L/10)).^2);
% c02 = @(x) (exp(100 * sin(x - 2)) - 1) ./ (exp(100 * sin(x - 2)) + 1);
c02 = @(x) (exp(100 * sin(x - 2)) - 1) ./ (exp(100 * sin(x - 2)) + 1);


kk = [0:N/2,-N/2+1:-1]';
omega = (2 * pi / L);
omegak = omega .* kk;



            % Fourier coefficient and Fourier series for c1 and c2

c01hat = fft(c01(xn));
c02hat = fft(c02(xn));

fc01 = ifft(c01hat);
fc02 = ifft(c02hat);
           % Plot the initial condition and hold it on the figure


%%
            % Part a: Pure advection (nu = 0, gamma = 0)
            % for u = 0.2, 0.5, and -0.5
            % use s = 0, and initial value c01

figure(1)
plot(xn, fc01, 'b')
title('Pure Advection')
hold on
uu = [0.2, 0.5, -0.5]';

for i = 1 : 3
    u = uu(i);
    a1 = -u;

    alphak = a1 * j * omegak;
    chat = c01hat .* exp(alphak * T);
    Cadvection = ifft(chat);

    plot(xn, Cadvection,'o--')
    hold on
end
legend('c01', 'u = 0.2' , 'u = 0.5', 'u = -0.5')



%%
            % Part b: Pure diffusion (u = 0, gamma = 0)
            % for nu = 0.01, 0.1, -0.01
            % use s = 0, and initial value c01


figure(2)
plot(xn, fc01, 'b')
hold on


vv = [0.01, 0.1, -0.01];
for i = 1 : 2 %for v = -0.01 the code is not working
    v = vv(i);
    a2 = v;

    alphak = a2 * (j * omegak) .^ 2;
    chat = c01hat .* exp(alphak * T);
    Cdiffusion = ifft(chat);

    plot(xn, Cdiffusion, 'o--')
    hold on

end
    

%%
            % Pure Advection
            % u = 0.5, s = 0,
            % Initial datum =  C02
            % Source value = 0

u = 0.5;
a1 = -u;

alphak = a1 * (omegak .* j);

chat = c02hat .* exp(alphak .* T);
Cadvection2 = ifft(chat);

figure(3)
plot(xn, fc02 ,'b', xn, Cadvection2, 'r--')
title('Pure Advection part C')
legend('co2', 'u = 0.5')

%%
            % Pure diffusion
            % v = 0.05, s = 0
            % Initial datum C02
            % Source value = 0

v = 0.01;
a2 = v;

alphak = a2 * (j * omegak) .^ 2;
chat = c02hat .* exp(alphak .* T);
Cdiffusion2 = ifft(chat);

figure(4)
plot(xn, fc02,'b', xn, Cdiffusion2,'r--')

%% 
            % Part e
            % Pure advection
            % u = 0.5, -0.5
            % Initial datum = 0
            % Source value S1

s1 = @(x) exp( -( (x - L/3) / (L/20)) .^ 2);
s1hat = fft(s1(xn));
fs1 = ifft(s1hat);

figure(5)
plot(xn, fs1, 'k')
hold on
uu = [0.5, -0.5]';

for i = 1 : 2
    u = uu(i);
    a1 = -u;

    alphak = a1 * (j * omegak);
    chat = zeros(N, 1);
    chat(1) = s1hat(1) * T;
    chat (2:N) = (exp(alphak(2:N) * T) -1) .* s1hat(2:N)./ alphak(2:N);
    cadvection3 = ifft(chat);

    plot(xn, cadvection3, 'd--')
    hold on
end
title('Pure advection with Source datum = 0')
legend('Source s1', 'u = 0.5', 'u = -0.5')

%% 
            % Part f
            % Pure diffusion
            % v = 0.5
            % Initial datum = 0
            % Source value S1

s2 = @(x) exp( -( (x - L/2) / (L/20)) .^ 2);
s2hat = fft(s2(xn));
fs2 = ifft(s2hat);

figure(6)
plot(xn, fs2, 'k')
hold on

v = 0.05;
a2 = v;

alphak = a2 * (j * omegak) .^ 2;
chat = zeros(N, 1);

chat (1) = s2hat(1) * T;
chat (2 : N) = s2hat(2 : N) .* (exp(alphak(2:N) * T) -1) ./ alphak(2 : N);
cdiffusion3 = ifft(chat);

plot(xn, cdiffusion3, 'ro-')
title('Pure diffusion with Initial datum = 0')
legend('source 2', 'ni = 0.5')

%%
            % Part g
            % advection diffusion
            % v = 0.5
            % Initial datum = c1
            % Source value S1
            % gamma = -0.01


s1 = @(x) exp( -( (x - L/3) / (L/20)) .^ 2);
s1hat = fft(s1(xn));
fs1 = ifft(s1hat);

u = 0.2;
v = 0.01;
gamma = - 0.01;

a0 = -gamma;
a1 = -u;
a2 = v;

alphak = a0 + a1 * j * omegak + a2 *(j * omegak).^ 2;

chat(1) = c01hat(1) + s1hat(1) * T;

chat(2:N) = c01hat (2:N) .* exp(alphak(2:N)*T) ...
   + s1hat(2:N) .* (exp(alphak(2:N) * T)-1) ./ alphak(2:N);

caddiff = ifft(chat);

figure(7)
plot(xn,fc01, 'b', xn, fs1, 'k', xn, caddiff, '*--')
legend('c02', 'source 1', 'CadvectionDiffusion')
title('adv & diff')

