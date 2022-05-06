%if we have periodic boundary condition



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



c0per = c0(xper);
c0hat = fft(c0per);

omega = 2 * pi / L;
kk = [0 : N/2, -N/2 + 1 : -1]';
alphak = ni * (j * omega * kk) .^ 2;
chat = exp(alphak * T) .* c0hat;

cper= ifft(chat)


dt * ni / (dx.^2)


% in case of a periodic boundary condition we can run this on the mesh that
% was for boundary condition, (the method that fourier method is using
cold = c0(xper);
cnew = cold;
t = 0;


% the first and the last points in the loop would be slightly different in
% this case. because the last point is equal to the first point

for n = 1 : nstep

    cnew(1) = cold(1) + dt * ni * (cold(2) - 2 * cold(1) + cold(N))/dx^2; %+dt*stime(xper(1), t))
    % here instead of g0(t) we put cold(N) because in the first point the
    % previous point we used to get the first point is exactly the last
    % point of cold

    for i = 2 : N-1 %also this goes one step more in periodic 
        
        cnew(i) = cold(i) + dt * ni * (cold(i+1) - 2 * cold(i) + cold(i-1))/(dx^2);%+dt*stime(xper(i), t))
        
    end
    
    %also first term here changes to cnew (N) it used to be cnew(N-1)
    cnew(N) = cold(N) + dt * ni * (cold(1) - 2 * cold(N) + cold(N-1))/(dx^2); %+dt*stime(xper(i), t))
    
    %here the last point in the periodic boundary condition is one point
    %after that of a non periodic boundary condition, so instead of N-1s we
    %have Ns now, and instead of gL(t) we have cold(1)
    cold = cnew;
    t = t + dt;
   
end

plot(xper, cnew, 'ro', xper, c0(xper), 'b', xper, cper(1:N), 'kd')

err = abs(cnew - cper(1:N));
figure()
plot(xper, err, 'bo')
norm(err, 2)/ norm(xper(1:N), 2)
norm(err, inf) / norm(cper(1:N), inf)