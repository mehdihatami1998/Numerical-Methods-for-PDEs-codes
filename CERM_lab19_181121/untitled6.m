           % Pure Advection
            % u = 0.5, s = 0,
            % Initial value  C02

clc
clear all
close all
format long

L = 10;
T = 5;
N = 200;

h = L / N;
x = [0 : h : L-h]';

c02 = @(x) (exp(100 * sin(x - 2)) - 1) ./ (exp(100 * sin(x - 2)) + 1);

c02hat = fft(c02(x));
fc02 = ifft(c02hat);

kk = [0 : N/2, (-N / 2 + 1) : -1]';
omega = 2 * pi / L;
omegak = kk * omega;

u = 0.5;
a1 = -u;

alphak = a1 * (omegak * j);

chat = c02hat .* exp(alphak * T);
Cadvection2 = ifft(chat);

plot(x, fc02,'b', x, Cadvection2, 'r--')
