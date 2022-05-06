% CERM_Lab_151121
% Exercise 3


clc
clear all
close all
format long

            % Defining the initial conditions
N = 100;
T = 20;
h = T / N;

x = [0 : h : T]';



for nmode = 1 : 200
    y= exp(-abs(x-5) .^ 2.5);
    yhat = fft(y) / N;
    figure(1)
    plot(x, abs(yhat).^2, 'ro');

    mask = ones(N, 1);
    mask(nmode + 2 : N - nmode) = 0;

    spectrum = abs( yhat(1 : N) ) .^ 2;

    numerator = sum(spectrum .* mask);
    denominator = sum(spectrum);
    ans = numerator / denominator;

    if ans >= 0.97
        break
    end

    i = i + 1;
end

nmode
ans