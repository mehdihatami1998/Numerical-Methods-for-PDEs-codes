% CERM_lab_21121
% Exercise 2
clear all
close all
clc



%%%%%%%%%%%%%%%%%%%% when diffusion coefficient is large enough, this
%%%%%%%%%%%%%%%%%%%% method would face efficienty issues. (e.g., when we
%%%%%%%%%%%%%%%%%%%% make the diffusity coefficient 3-5 times larger *3),
%%%%%%%%%%%%%%%%%%%% and how increasing the time steps can solve that
%%%%%%%%%%%%%%%%%%%% problem?,

% exact solution
yex = @(x, t) exp(-t) .* sin(8*x);

L = 2*pi;
T = 1;


N = 50;
dx = L / N;

M = 4000;
dt = T / M;


% source function
s = @(x, t) -sin(8*x) * exp(-t) - 8 .* (cos(8*x) .* exp(-t)./(t+1)) ...
            + 64 .* (sin(8.*x) .* exp(-t) .* (x+1)/(t+1));


ni = @(x, t) (1+x) / (1+t);

% Diriclet boundary condition

g0 = @(t) 0;
gL = @(t) 0;


% Initial condition 
c0 = @(x) sin(8*x);


% internal nodes for dirichlet boundary condition
% xhalf for nuvar (not constant nu)

xin = [dx : dx : L-dx]';
xhalf = [dx * 0.5 : dx : L - dx * 0.5]';

          % check stability

% if dt * v / (dx.^2) < 0.5
%     display ('Stable Condition')
% else
%     display ('Unstable Conditon')
% end

            % stability check charlie's code



t = 0;
cold = c0 (xin);
cnew = cold;

for n = 1:M

    cnew(1) = cold(1) + (dt/dx.^2) * (ni(xhalf(2), t) * (cold(2) - cold(1))...
                                    - ni(xhalf(1), t) * (cold(1) - g0(t))) ...
                                    + dt * s(xin(1), t);

    for i = 2 : N-2

        cnew(i) = cold(i) + (dt/dx.^2) * (ni(xhalf(i+1), t)*(cold(i+1) - cold(i))...
                                        - ni(xhalf(i), t)* (cold(i) - cold(i-1))) ...
                                        + dt * s(xin(i), t);

    end

    cnew(N-1) = cold(N-1) + (dt/dx.^2) * (ni(xhalf(N), t)*(gL(t) - cold(N-1))...
                                    - ni(xhalf(N-1), t)* (cold(N-1) - cold(N-2))) ...
                                    + dt * s(xin(N-1), t);
    t = t+dt;
    cold = cnew;

end
figure(1)
plot(xin, cnew, 'r*--',xin, c0(xin), 'b')
legend('after time T', 'initial condition')

figure(2)
plot(xin, cnew, 'r-*');
hold on

error1 = abs(yex(xin, T) - cnew);
norm2error1= norm(error1,2) * dx;

norminferror1 = norm(error1, inf)
%%%%%% second solution


% CERM_lab_21121
% Exercise 2


% exact solution
yex = @(x, t) exp(-t) .* sin(8*x);

L = 2*pi;
T = 1;


N = 50;
dx = L / N;

M = 4000;
dt = T / M;


% source function
s = @(x, t) -sin(8*x) * exp(-t) - 8 .* (cos(8*x) .* exp(-t)./(t+1)) ...
            + 64 .* (sin(8.*x) .* exp(-t) .* (x+1)/(t+1));


ni = @(x, t) (1+x) / (1+t);

% Diriclet boundary condition

g0 = @(t) 0;
gL = @(t) 0;


% Initial condition 
c0 = @(x) sin(8*x);


% first solution
N = 100;
dx = L /N;

M = 8000;
dt = T / M;


% internal nodes for dirichlet boundary condition
% xhalf for nuvar (not constant nu)

xin = [dx : dx : L-dx]';
xhalf = [dx * 0.5 : dx : L - dx * 0.5]';

          % check stability

% if dt * v / (dx.^2) < 0.5
%     display ('Stable Condition')
% else
%     display ('Unstable Conditon')
% end

            % stability check charlie's code



t = 0;
cold = c0 (xin);
cnew = cold;

for n = 1:M

    cnew(1) = cold(1) + (dt/dx.^2) * (ni(xhalf(2), t) * (cold(2) - cold(1))...
                                    - ni(xhalf(1), t) * (cold(1) - g0(t))) ...
                                    + dt * s(xin(1), t);

    for i = 2 : N-2

        cnew(i) = cold(i) + (dt/dx.^2) * (ni(xhalf(i+1), t)*(cold(i+1) - cold(i))...
                                        - ni(xhalf(i), t)* (cold(i) - cold(i-1))) ...
                                        + dt * s(xin(i), t);

    end

    cnew(N-1) = cold(N-1) + (dt/dx.^2) * (ni(xhalf(N), t)*(gL(t) - cold(N-1))...
                                    - ni(xhalf(N-1), t)* (cold(N-1) - cold(N-2))) ...
                                    + dt * s(xin(N-1), t);
    t = t+dt;
    cold = cnew;

end


plot(xin, yex(xin, T), 'k', xin, cnew, 'b--o')
legend('Euler M = 4000', 'Exact Solution', 'Euler = 8000')


error2 = abs(yex(xin, T) - cnew);
norm2error2= norm(error2,2) * dx;
norminferror2 = norm(error2, inf)

p_emp = -log2(norm2error2/norm2error1)
p_emp1 = -log2(norminferror2/norminferror1)
