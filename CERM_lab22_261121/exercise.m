% separation of variables



clear all
clc
close all
format long



L = 10;
T = 3;
N = 100;
dx = L / N;

% intervals
xper = [0 : dx : L - dx]';
xin = [dx : dx : L - dx]';

c0 = @(x) 10 * exp(-((x-L/2)/(L/10)).^2);
c0hat = fft(c0(xper));
fc0 = ifft(c0hat);

kk = [0 : N/2, -N/2 + 1 :-1]';
omega = 2 * pi / L;
omegak = omega * kk;

v = 0.1;
a2 = v;

alphak = a2 *(omegak * j).^2;

chat = c0hat .* exp(alphak .* T);
fc = ifft(chat);

figure()
plot(xper, fc0, 'r', xper, fc, 'k--d')




%%
% implicit euler method

% clc
% clear all
% close all 
% format long
% 

L = 10;
T = 3;
N = 100;
M = 50;
dx = L / N;
dt = T / M;


ni = 0.1;


cold = c0(xper);
cnew = cold;
dtdx2 = dt/dx.^2;

t = 0;
A = zeros(N-1);
b = zeros(N-1, 1);


for n = 1 : M

    b = cold + 0;

    b(1, 1) = b(1, 1) - ni * dtdx2;

    b(N-1, 1) = b(N-1, 1) - ni * dtdx2;

    A = diag( ones(N-1) + 2 * ni * dtdx2);
    
    A = A - diag( ones(N-1) * ni * dtdx2 ,1) - diag(ones(N-1) * ni * dtdx2, -1);

    cnew = A \ b;

    cold = cnew;
     
    t = t + dt;
end

figure(1)
plot(xin, cnew)

