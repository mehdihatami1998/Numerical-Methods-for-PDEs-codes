clear all
close all
format long
clc

L = 10;
T = 2;
ni = 0.01;

c0 = @(x) sin(2 * pi * x / L);


% boundary condition:
g0 = @(t)  0;
gL = @(t)  0;



N = 50;
nstep = 100;


dx = L / N;
dt = T / nstep;
% this is the mesh:
x = [0 : dx : L]';


% this is the internal nodes:
xin = [dx : dx : L - dx]'; % for dirichlet boundary condition
xper = [0 : dx : L - dx]'; % for periodic boundary condition




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% this part added if periodic
c0per = c0(xper);
c0hat = fft(c0per);

omega = 2 * pi / L;
kk = [0 : N/2, -N/2 + 1 : -1]';
alphak = ni * (j * omega * kk) .^ 2;
chat = exp(alphak * T) .* c0hat;

cper= ifft(chat)

% this is computing a reference solution, and it works only if our function
% is periodic, in that case it gives accurate results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check stability restrictions:
% this number has to be smaller than 0.5
dt * ni / (dx.^2)

cold = c0(xin);
% so the idea is that we need to have a double loop,
% one loop for the time and one loop for the space.


% we put cnew = cold to get the size of cold for the cnew
cnew = cold;
t = 0;
for n = 1 : nstep

    cnew(1) = cold(1) + dt * ni * (cold(2) - 2 * cold(1) + g0(t))/dx^2;
    % above we used the same formula as we have in the line below, the
    % difference is between the term cold(i-1) that here we have a boundary
    % condition which we have to implement at time 

    % also the same thing happens for the last node that I'll write after
    % the for iteration for cn(N-1)

    for i = 2 : N-2 
        
        cnew(i) = cold(i) + dt * ni * (cold(i+1) - 2 * cold(i) + cold(i-1))/(dx^2);
        %above is the implementation of formula 11.11 in case of a constant
        %ni
    end
    
    cnew(N-1) = cold(N-1) + dt * ni * (gL(t) - 2 * cold(N-1) + cold(N-2))/(dx^2);
    % here unlike the cnew(1) in which we should use boundary condition in
    % the first point, here as the last point, we should change
    % cold(i+1)term to the boundary condition term
    cold = cnew;
    t = t + dt;
   
end

plot(xin, cnew, 'ro', xin, c0(xin), 'b', xin, cper(2:N), 'kd')

err = abs(cnew - cper(2:N));
figure()
plot(xin, err, 'bo')
norm(err, 2)/ norm(xper(2:N), 2)
norm(err, inf) / norm(cper(2:N), inf)