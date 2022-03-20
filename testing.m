%% %Newton method
close all
clear all
format long
f=@(x) x.^3-x.^2-2.*x;
df=@(x) 3*x.^2 - 2.*x -2;
x=[-3:0.01:3]';
y=f(x);
plot(x,y)
grid on 
%answers = -1, 0, 2
x0=-0.8;
tolerance=10.^-7;
maxIter=200;
for i=1:maxIter
    i
    delta=-f(x0)/df(x0);
    if abs(f(x0))<tolerance
        break
    end
    x0=x0+delta;
end
x0
%% %fixed point method
close all
clear all
format long
f=@(x) x.^3-x.^2-2.*x;
df=@(x) 3*x.^2 - 2.*x -2;
phi=@(x) 0.5*(x.^3-x.^2);
dphi=@(x) 0.5.*(3.*x.^2-2.*x); 
x=[-3:0.01:3]';
y=f(x);
plot(x,abs(dphi(x)), x,ones(size(x)),'k*')
grid on 
figure(2)
plot(x,y)

%
figure(3)
tiledlayout(2,1)
nexttile
plot(x,abs(dphi(x)), x,ones(size(x)),'k')
grid on
nexttile
plot(x,y)
grid on
%answers = -1, 0, 2
x0=1.6;
tolerance=10.^-7;
maxIter=200;
for i=1:maxIter
    i
    x1=phi(x0);
    if abs(f(x1))<tolerance
        break
    end
    x0=x1
end
x1