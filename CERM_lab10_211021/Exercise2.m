
%%
            % Solving a non linear equation
            % Secand Method
close all
clear all
format long

f = @(x) x.^ 5 - 8.* x.^ 4 + 21.* x.^ 3 - 14.* x.^ 2 -20.* x + 24;
df = @(x) ( - 3* x^ 2+ 5* x ) / (x - 1) ^2+ (6* x- 5) / (x - 1);

x = [1: 0.01: 3.5]';
y = f(x);

plot(x, y, x, zeros(size(x)), 'k*');
grid on 


            % in this plot we see there are 2 answers for this equation
            % answers are approximately around 2 and 3
            % we want to find the one which is around 3, and for using 
            % secant method we need to have 2 initial guesses, so we choose 
            % x1 = 2.7 and x2 = 3.15 as our initial guess to use in Newton's method

x2 = 2.7;
x1 = 3.15  ;

tolerance = 10e-11;
maxIter=200;


for i = 1 : maxIter
    delta = -f(x2).* (x1-x2) / ( f(x1) - f(x2) );

    if abs(delta) < tolerance
        break
    end

   x1 = x2;
   x2 = x2+ delta;
end

x2
sprintf (" the number of iterations for Newton's method = %.15g ", i)

%%
            % Newton method

close all
clear all
format long

f = @(x) x.^ 5 - 8.* x.^ 4 + 21.* x.^ 3 - 14.* x.^ 2 -20.* x + 24;
df = @(x) ( - 3* x^ 2+ 5* x ) / (x - 1) ^2+ (6* x- 5) / (x - 1);

x = [1: 0.01: 3.5]';
y = f(x);

plot(x, y, x, zeros(size(x)), 'k*');
grid on 

            
            % in this plot we see there are 2 answers for this equation
            % answers are approximately around 2 and 3
            % we want to find the one which is around 3, and for using 
            % Newton method we need to have only one initial guesses, so we choose 
            % x1 = 2.85 as our initial guess to use in Newton's method

x0=2.85;

tolerance = 10e-10;
maxIter=200;

for i = 1 : maxIter

    delta = -( f(x0)/ df(x0) );
    if abs( delta/(x0) )<tolerance
        break
    end
    x0 = x0+ delta;
end

x0

sprintf (" the number of iterations for Newton's method = %.15g ", i)
