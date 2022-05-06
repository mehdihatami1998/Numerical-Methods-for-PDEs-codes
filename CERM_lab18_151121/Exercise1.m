% CERM_Lab_151121
% Exercise 1

clc
clear all
close all
format long

            % Defining the parameters
N = 100;        % repeat the code with N = 10, 20, 40, 50
L = 10;
h = L / N;
x = [0 : h : L-h]';


            % Defining the functions

f1 = sin(4 * pi / L *x);

f2 = exp(-((x - L/2) / (L/10)).^2);

f3 = exp(- abs((x - L/2) / (L/20)));

f4 = exp(-((x - L/2) / (L/10)).^2) + 0.001 * randn(N, 1);

f5 = exp(-((x - L/2) / (L/10)).^2) + 0.1 * randn(N, 1);


 


            % Defining the fourier coefficients

fhat1 = fft(f1)/N;
fhat2 = fft(f2)/N;
fhat3 = fft(f3)/N;
fhat4 = fft(f4)/N;
fhat5 = fft(f5)/N;

            % Plot spectrum for each function

figure(1)
plot([1:N/2], abs(fhat1(1:N/2)).^2, 'ro')

figure(2)
plot([1:N/2], abs(fhat2(1:N/2)).^2, 'ro')

figure(3)
plot([1:N/2], abs(fhat3(1:N/2)).^2, 'ro')


figure(4)
plot([1:N/2], abs(fhat4(1:N/2)).^2, 'ro')

figure(5)
plot([1:N/2], abs(fhat5(1:N/2)).^2, 'ro')


