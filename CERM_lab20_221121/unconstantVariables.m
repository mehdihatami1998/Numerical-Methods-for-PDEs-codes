%if we have periodic boundary condition



clear all
close all
format long
clc

L = 10;
T = 4;
ni = 0.1;

c0 = @(x) sin(2 * pi * x / L);


% boundary condition:
g0 = @(t)  0;
gL = @(t)  0;


nivar = @(x, t) ni * (ones(size(x)) + cos(2 * pi * x / L).^2) + (1+ 0.2 * sin(t));
% based on the formula in the lecture notes, this "nivar" needs to be
% computed at the nodes at Xhalf, so we define the xhalf period as well.
% these points are centers of the mesh


N = 50;
nstep = 100;


dx = L / N;
dt = T / nstep;
% this is the mesh:
x = [0 : dx : L]';
xhalf = [dx * 0.5 : dx : L - dx * 0.5]';


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
cold = c0(xin);
cnew = cold;
t = 0;


% the first and the last points in the loop would be slightly different in
% this case. because the last point is equal to the first point

for n = 1 : nstep

    cnew(1) = cold(1) + dt * ( nivar(xhalf(2) , t) * (cold(2) -  cold(1)) ...
                                 - nivar(xhalf( 1)    , t) * (cold(1)   - g0(t))) /(dx^2); %+dt*stime(xper(1), t))
    

    for i = 2 : N-2 %for periodic "N-1" for dirichlet "N-2"
        
        cnew(i) = cold(i) + dt * ( nivar(xhalf(i + 1) , t) * (cold(i+1) -  cold(i)) ...
                                 - nivar(xhalf( i)    , t) * (cold(i)   - cold(i-1))) /(dx^2); %+dt*stime(xper(i), t))
         
    end
    
    %also first term here changes to cnew (N) it used to be cnew(N-1)
    cnew(N-1) = cold(N-1) + dt * ( nivar(xhalf(N) , t) * (gL(t) -  cold(N-1)) ...
                                 - nivar(xhalf( N-1), t) * (cold(N-1)   - cold(N-2))) /(dx^2); %+dt*stime(xper(i), t))
    
    %here the last point in the periodic boundary condition is one point
    %after that of a non periodic boundary condition, so instead of N-1s we
    %have Ns now, and instead of gL(t) we have cold(1)
    cold = cnew;
    t = t + dt;
   
end

plot(xin, cnew, 'r--', xin, c0(xin), 'b', xin, cper(2:N), 'kd')

% err = abs(cnew - cper(1:N));
% figure()
% plot(xper, err, 'bo')
% norm(err, 2)/ norm(xper(1:N), 2)
% norm(err, inf) / norm(cper(1:N), inf)