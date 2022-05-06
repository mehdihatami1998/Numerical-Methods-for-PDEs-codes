clear all
close all
clc 
format long


L = 10;
N = 400;
T = 1;

u = 4

a1 = -u; % a1*T it says how much and in which direction we should go to take 
        % the sampelling point. if it is negative as in this case, it will
        % go forward, and if it's positive, it will go in negative
        % direction

         
        % a1 is the opposite of velocity, and in many exercises we will be
        % asked to interpret the opposite of velocity which is u

nmodes = 100;

h = L / N;

x = [0 : h : L-h]';

c0 = @(x) exp(-(x - L * 0.5) .^2 / (L * 0.1) .^2);


c0hat = fft(c0(x));

% mask = ones(N,1);
% mask(nmodes+2:N-nmodes) = 0
% 
%advection equation

omega = 2 * pi / L;
%  kk = [0;[1:N/2]';[-(N/2)+1):-1]'];

kk = [0 : N/2,(-N/2+1):-1]';
omegak = omega * kk ; % j is the imaginary unit


% dealing with zero source K
alphak =  a1 * omegak * j; % this is the order that matlab makes fft


chat = c0hat .* exp(alphak*T);



% solution of time t
% c = ifft(chat.*mask);

c = ifft(chat);

plot(x, c0(x), 'b', x, c, 'r')

