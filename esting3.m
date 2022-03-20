%%
%Newton method
close all
clear all 
format long
f=@(x) cos(5*x).*x.^(3/7);
df=@(x) (3*cos(5.*x))./(7.*x.^(4/7)) - 5*x.^(3/7).*sin(5.*x);
x=[0:0.01:2]';
y=f(x);
plot(x,y) %0.3, 0.9, 1.65
grid on
x=1.6;
tolerance=10e-11;
maxIter=200;
for i=1:maxIter
    i
    delta=-f(x)/df(x);
    if abs(delta)<tolerance
        break
    end
    x=x+delta
end
x %0.314159 , 0.942477, 1.5707963
%%
%bisection method
close all
clear all 
format long
f=@(x) cos(5*x).*x.^(3/7);
df=@(x) (3*cos(5.*x))./(7.*x.^(4/7)) - 5*x.^(3/7).*sin(5.*x);
x=[0:0.01:2]';
y=f(x);
plot(x,y) %0.3, 0.9, 1.65
grid on
a=0.7;
b=1;
xm=(a+b)/2;
tolerance=10e-11;
maxIter=200;
for i=1:maxIter
    i
    if f(a).*f(xm)<0
        b=xm;
    else
        a=xm;
    end
    xm=(a+b)./2;
    err=(a-b)/2
        if abs(err) < tolerance
            break
        end
end
xm 

%Newton method finishes operating at 4th operation to get the answer with a
%tolerance of 10e-11, bisection method takes 34 iteration for the same
%operation.