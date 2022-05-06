% writing  the for loop in vector form:
%%%%%%%%%%%%%% comment:
% So when we have a large number of timesteps, the explicite method can 
% work as good as the implicit method, say 500 timesteps, and there won't
% be any significant difference in their errors, but as soon as we decrease
% the number of timesteps, (to a degree that the explicit method is still
% stable) we would see that the error of explicit method will increase
% significantly, and if we continue decreasing the number of timesteps, the
% stability condition for the explicit method would not be satisfied but
% the implicit method would continue working even as small number for time
% steps as 5 or 10.


% this is with periodic boundary condition
% the difference is that in periodic boundary condition the number of
% unknowns is not N-1, it is N. so we have N unknowns in periodic B.C.
% mode. because not only do we have the internal points, we also have one
% of the boundary points 



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


N = 200;
nstep = 500;


dx = L / N;
dt = T / nstep;
% this is the mesh:
x = [0 : dx : L]';
xhalf = [dx * 0.5 : dx : L - dx * 0.5]';


xin = [dx : dx : L - dx]'; % for dirichlet boundary condition
xper = [0 : dx : L - dx]'; % for periodic boundary condition



c0per = c0(xper);
c0hat = fft(c0per);
shat = fft(s(xper));

omega = 2 * pi / L;
kk = [0 : N/2, -N/2 + 1 : -1]';
alphak = ni * (j * omega * kk) .^ 2;
chat = zeros(N,1);

chat(1,1) = c0hat(1, 1) + T * shat(1, 1);
chat(2:N, 1) = exp(alphak(2:N, 1) * T) .* c0hat(2:N, 1)...
            + shat(2:N, 1) .* (exp(alphak(2:N, 1) * T)-1) ./ alphak(2:N, 1);

cper= ifft(chat);


dt * ni / (dx.^2);

cold = c0(xper);
cnew = cold;
cnew1 = cold; 
cold1 = cold;
t = 0;

A = zeros(N-1);
% b = zeros(N-1, 1);
dtdx2=dt/dx^2;


for n = 1 : nstep
    
    b = cold1 + s(xper) * dt;


    % the following two terms are only needed in the case of explicit
    % boundary condition, so we don't need them here with periodic boundary
%     % condition. so I commented them out.
%     b(1,1) = b(1, 1) - dtdx2 * nivar(xhalf(1),t + dt) * g0(t+dt);
%     b(N-1, 1) = b(N-1, 1) - dtdx2 * nivar(xhalf(N), t+dt) * gL(t+dt);
    


    A = diag(ones(N-1, 1) + nivar(xhalf(2:N), t + dt) * dtdx2 + nivar(xhalf(1:N-1), t+dt) * dtdx2)...
        - diag(ones(N-2,1) .* nivar(xhalf(3:N), t+dt) * dtdx2, 1)...
        - diag(ones(N-2,1) .* nivar(xhalf(1:N-2), t +dt) * dtdx2, -1); 

    cnew1 = A\b;

    cold1=cnew1;

    t = t + dt;
   
end

plot(xin, cnew1, 'r*', xin, cper(2:N), 'kd')
legend('cnew', 'cper(2:N)')

% figure()
% plot(xin, cper(2:N))
err = abs(cnew1 - cper(2:N));

norm(err, 2) * sqrt(dx)/ (norm(xper(1:N), 2) *sqrt(dx))
norm(err, inf) / norm(cper(1:N), inf)