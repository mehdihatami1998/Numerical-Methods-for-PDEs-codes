% writing  the for loop in vector form:


%if we have periodic boundary condition



clear all
close all
format long
clc

L = 10;
T = 6;
ni = 0.1;

c0 = @(x) sin(2 * pi * x / L);


% boundary condition ( homogenous )
g0 = @(t)  0;
gL = @(t)  0;


nivar = @(x, t) ni * (ones(size(x)) );%+ cos(2 * pi * x / L).^2) + (1+ 0.2 * sin(t));

s = @(x) exp(-(x-L*0.3).^2/(L*0.05)^2);


N = 50;
nstep = 100;


dx = L / N;
dt = T / nstep;
% this is the mesh:
x = [0 : dx : L]';
xhalf = [dx * 0.5 : dx : L - dx * 0.5]';


xin = [dx : dx : L - dx]'; % for dirichlet boundary condition
xper = [0 : dx : L - dx]'; % for periodic boundary condition



c0per = c0(xper);
c0hat = fft(c0per);

omega = 2 * pi / L;
kk = [0 : N/2, -N/2 + 1 : -1]';
alphak = ni * (j * omega * kk) .^ 2;
chat = exp(alphak * T) .* c0hat;

cper= ifft(chat);


dt * ni / (dx.^2);

cold = c0(xin);
cnew = cold;
cnew1 = cold; 
cold1 = cold;
t = 0;

A = zeros(N-1);
b = zeros(N-1, 1);
dtdx2=dt/dx^2;


for n = 1 : nstep
    
    b = s(xin) * dt;
    b(1,1) = b(1, 1) + dtdx2 * nivar(xhalf(1),t) * g0(t);

    b(N-1, 1) = b(N-1, 1) + dtdx2 * nivar(xhalf(N), t) * gL(t);


    A = diag(ones(N-1, 1) - nivar(xhalf(2:N), t) * dtdx2 - nivar(xhalf(1:N-1), t) * dtdx2);
    A = A + diag(nivar(xhalf(3:N), t) * dtdx2, 1) + diag(nivar(xhalf(1:N-2), t) * dtdx2, -1);


%     cnew(1) = cold(1) + dt * ( nivar(xhalf(2) , t) * (cold(2) -  cold(1)) ...
%                                  - nivar(xhalf( 1) , t) * (cold(1)   - g0(t))) /(dx^2) + dt * s(xper(1));


    cnew(1) = cold(1) * (1 - nivar(xhalf(2), t) * dtdx2 - nivar(xhalf(1), t) * dtdx2)...
            + cold(2) * nivar(xhalf(2), t) * dtdx2 ...
            +g0(t) * nivar(xhalf(1), t) * dtdx2 + s(xper(1))*dt;
            + cold()


    for i = 2 : N-2 
        
%         cnew(i) = cold(i) + dt * ( nivar(xhalf(i + 1) , t) * (cold(i+1) -  cold(i)) ...
%                                  - nivar(xhalf( i)    , t) * (cold(i)   - cold(i-1))) /(dx^2)+ dt* s(xper(i));
         
        cnew(i) = cold(i) * (1 - nivar(xhalf(i + 1), t) * dtdx2 - nivar(xhalf(i), t)*dtdx2)...
                 + cold(i+1) * (nivar(xhalf(i+1), t) * dtdx2)...
                 + cold(i-1) * (nivar(xhalf(i), t) * dtdx2) + dt * s(xper(i))


                    % this inner loop here can be written as the
                    % multipllication of a vecot plus a constant term in
                    % order to get things easier in systems, we need:
                    % Cnew = A * cold + b instead of the whole loop here

    end
    
    cnew(N-1) = cold(N-1) + dt * ( nivar(xhalf(N) , t) * (gL(t) -  cold(N-1)) ...
                                 - nivar(xhalf( N-1), t) * (cold(N-1)   - cold(N-2))) /(dx^2) + dt*s(xper(i));
    
    cnew1 = A * cold1 + b;
    cold = cnew;
    cold1 = cnew1;
    t = t + dt;
   
end

plot(xin, cnew, 'ro', xin, cnew1, 'b*', xin, cper(2:N), 'kd')
legend('cnew', 'c0(xin)', 'cper(2:N)')

% err = abs(cnew - cper(1:N));
% figure()
% plot(xper, err, 'bo')
% norm(err, 2)/ norm(xper(1:N), 2)
% norm(err, inf) / norm(cper(1:N), inf)