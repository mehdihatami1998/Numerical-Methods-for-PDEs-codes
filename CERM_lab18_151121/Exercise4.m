% CERM_Lab_151121
% Exercise 3


clc
clear all
close all
format long


            % Defining the function and the inital condition
L = 8;
N = 160;
h = L / N;

x = [0 : h : L-h]';
y1 = @(x) (x/L) .*sin((10 * pi .* x) / L);

y = y1(x);
dy1 = @(x) (1/L) .* sin((10 * pi .* x) / L) + ((10 * pi .* x) ./ (L .^ 2)) ....
     .* cos((10 * pi * x) / L);

dy = dy1(x);


%%%%%%%%%%%%%%%% Comments %%%%%%%%%%%%%%%%
% If we use fourier series to approximate a regular function, the error
% will go to zero extremely fast, so fast that there won't be a constant
% number for the value of empirical estimate of the convergence order, but
% after reaching the accuracy of machine epsilone, the P_emp will start
% decreasing again, because the error is so small that the truncation error
% will play a role in that calculation.

%%%%%%%%%%%%%%%% Comments %%%%%%%%%%%%%%%%
% If a function is not periodic it's approximation would be like the
% approximation of a periodic function, the fourier series would
% approximate that as if it was periodic, with a discontinuoity at zero 


            % Defining masks

nmodes1 = 30;
nmodes2 = 60;

mask1 = ones(N, 1);
mask1(nmodes1 + 2 : N - nmodes1) = 0;

mask2 = (ones(N, 1));
mask2(nmodes2 + 2 : N - nmodes2) = 0;


figure(1)
plot(x,y,'r',x,dy,'b')


            % Fourier series  

yhat1=fft(y);
yy1 = ifft(yhat1 .* mask1);
yy2 = ifft(yhat1 .* mask2);


figure(2)
semilogy(x(1:N/2), (abs(yhat1(1:N/2))/N).^2,'ro')



            % Norms calculation

norm_inf_1 = norm(abs(y - yy1), inf);
norm_two_1 = norm(abs(y - yy1),  2)

norm_inf_2 = norm(abs(y - yy2), inf);
norm_two_2 = norm(abs(y - yy2),  2);

p_emp = -log2(norm_inf_2 / norm_inf_1)