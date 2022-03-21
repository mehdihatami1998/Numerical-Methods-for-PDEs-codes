%%
            % Solving a non linear equation
            % Newthon's Method
close all
clear all 
format long

f = @(x) cos(5* x).* x.^ (3/7);
df = @(x) ( 3* cos(5.* x) ) ./ ( 7.* x.^ (4/7) ) - 5* x.^ (3/7).* sin(5.* x);

            % making a mesh on x axis with the step of 0.01
x = [0: 0.01: 2]';
y = f(x);


plot(x, y, x, zeros(size(x)),'k*') 
grid on
            % in this plot we see there are 3 answers for this equation
            % answers are approximately around %0.3, 0.9, 1.65
            % we want to find the one which is around 1.6, so we choose 
            % x = 1.6 as our initial guess to use in Newton's method


x = 1.6;
tolerance = 10e-12;
maxIter = 200;

for i = 1:maxIter
    i;
    delta = -f(x) / df(x);
    if abs(delta) < tolerance
        break
    end
    x = x + delta;
end

sprintf (" the number of iterations for Newton's method = %.15g ", i)
x




%%
            % Bisection Method

close all
clear all 
format long

f = @(x) cos(5* x).* x.^ (3/7);
df = @(x) ( 3* cos(5.* x) )./ ( 7.* x.^ (4/7) ) - 5 * x.^ (3/7).* sin(5.* x);

               % making a mesh on x axis with the step of 0.01
x = [ 0: 0.01: 2]';
y = f(x);


plot(x, y, x, zeros(size(x)),'k*') 
grid on


            % in this plot we see there are 3 answers for this equation
            % answers are approximately around %0.3, 0.9, 1.65
            % we want to find the one which is around 1.6, so we choose 
            % x = 0.9 as our initial guess to use in Newton's method


a = 0.7;
b = 1;
xm = (a+b) / 2;

tolerance = 10e-12;
maxIter = 200;

for i = 1 : maxIter
    i;
    if f(a).* f(xm) < 0
        b = xm;
    else
        a = xm;
    end

    xm = (a+b)./ 2;
    err = (a-b)/2

        if abs(err) < tolerance
            break
        end
end
sprintf(" the number of iterations for Bisection method = %.15g ", i)
xm 

%Newton method finishes operating at 4th operation to get the answer with a
%tolerance of 10e-11, bisection method takes 31 iteration for the same
%operation.