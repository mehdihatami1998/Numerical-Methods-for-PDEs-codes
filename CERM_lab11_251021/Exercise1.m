            % Numerical Solutions of ODEs
            % Forward Euler Method
clc
clear all
close all
format long

f = @(t,y) -3* t.^2.* y;
yex = @(t) 3* exp(-t.^3);

y0 = 3;
T = 2;
h = 0.01;
N = T/ h;

tn = zeros(N+1,1);
un = zeros(N+1,1);
un(1,1) = y0;

for k = 1:N
    un(k+1,1) = un(k,1) + h* f( tn(k,1),un(k,1) );
    tn(k+1,1) = tn(k,1) + h;
end




            % Plotting Euler Method approximation vs the Exact answer
tmesh = [0: 0.01: T]';
figure(1)
plot(tmesh, yex(tmesh), 'r', tn, un, 'kd')
title("EULER Method vs. Exact function")

            % Absolute and Relatve Errors on the Interval

abs_err_EULER = abs( yex(tn) - un );
rel_err_EULER = abs_err_EULER/ yex(tn);

figure(2)
semilogy(tn,rel_err_EULER)
title("EULER Method vs. Relative Error")

            % Absolute and Relative Errors at Final Time

abs_error_endpoint_EULER = abs( yex(tn(end)) - un(end) )
rel_error_endpoint_EULER = abs( yex(tn(end)) - un(end)) / abs( yex( tn(end) ) )


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_EULER = norm( abs_err_EULER, inf )
rel_err_Norminf_EULER = norm( abs_err_EULER, inf ) / norm( yex(tn),inf )


%%
            % Heun Method
% clear all
% close all
% format long

f = @(t,y) -3* t.^2.* y;
yex = @(t) 3* exp(-t.^3);

y0 = 3;
T = 2;
h = 0.01;
N = T/ h;

tn = zeros(N+1,1);
un = zeros(N+1,1);
un(1,1) = y0;

for k = 1:N
    tn(k+1, 1) = tn(k,1) + h;
    f1 = f( tn(k,1), un(k,1) );
    f2 = f( tn(k+1,1), un(k,1)+ h*f1);
    un(k+1,1) = un(k,1) + 0.5* h* (f1+f2);
end




            % Plotting HEUN Method approximation vs the Exact answer
tmesh = [0: 0.01: T]';
figure(3)
plot(tmesh, yex(tmesh), 'r', tn, un, 'kd')
title("HEUN Method vs. Exact function")


            % Absolute and Relatve Errors on the Interval

abs_err_HEUN = abs( yex(tn) - un );
rel_err_HEUN = abs_err_HEUN/ yex(tn);

figure(4)
semilogy(tn,rel_err_HEUN)
title("HEUN Method Relative Error")


            % Absolute and Relative Errors at Final Time

abs_error_endpoint_HEUN = abs( yex(tn(end)) - un(end) )
rel_error_endpoint_HEUN = abs( yex(tn(end)) - un(end)) / abs( yex( tn(end) ) )


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_HEUN = norm( abs_err_HEUN, inf )
rel_err_Norminf_HEUN = norm( abs_err_HEUN, inf ) / norm( yex(tn),inf )

%% part c: leapfrog method:
clear all
close all
format long
f=@(t,y) -3*t.^2.*y;
y0=3;
T=2;
yex=@(t) 3*exp(-t.^3);
h=0.01;
N=T/h;
tn=zeros(N+1,1);
un=zeros(N+1,1);
un(1,1)=y0;
tn(2,1)=0.01;
f1=f(tn(1),un(1));
f2=f(tn(2),un(1)+h.*f1);
un(2,1)=y0+ 0.5*h*(f1+f2);
for k=2:N
    tn(k+1,1)=tn(k,1)+h;
    un(k+1,1)=un(k-1,1) + 2*h*f(tn(k,1),un(k,1));
end

tt=[0:0.01:T]';
figure(1)
plot(tt,yex(tt),'r',tn,un,'ko')
abs_err=abs(yex(tn)-un)
rel_err=abs_err/yex(tn);
figure(2)
semilogy(tn,abs_err)
abs_err_inf=norm(abs_err,inf)
rel_err_inf=norm(abs_err,inf)/norm(yex(tn),inf)
