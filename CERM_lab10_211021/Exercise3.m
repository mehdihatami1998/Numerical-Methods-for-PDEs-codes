%%
            % Solving a non linear equation
            % Newton Method
close all
clear all 
format long

f = @(x) x.^ 3- x.^ 2- 2.* x;
df = @(x) 3* x.^ 2 - 2.* x -2 ;

            % making a mesh on x axis with the step of 0.01

x = [-3: 0.01: 3]';
y=f(x);

plot(x,y)
grid on 

            % in this plot we see there are 3 answers for this equation
            % answers are approximately around -1, 0, 2
            % we want to find the one which is around 1.6, so we choose 
            % x = -0.8 as our initial guess to use in Newton's method


x0 = 0.3;

tolerance = 10e-8;
maxIter=200;


for i=1:maxIter

    delta = -f(x0)/ df(x0);
    
    if abs( f(x0) ) < tolerance
        break
    end

    x0 = x0+ delta;
end

x0
sprintf (" the number of iterations for Newton's method = %.15g ", i)


%% 
            %Fixed Point Method
close all
clear all
format long

f = @(x) x.^3 - x.^2 - 2.*x;
df=@(x) 3* x.^2 - 2.*x -2;
phi = @(x) 0.5* (x.^3 - x.^2);
dphi = @(x) 0.5.* (3.* x.^2 - 2.*x); 


            % making a mesh on x axis with the step of 0.01 over the given
            % interval

x = [ -3: 0.01: 3]';
y = f(x);

            % Based on the theory, we know that fixed point method only can
            % provide us with the answers which are in the interval that
            % the first derivative of phi function is below 1



figure(1)
tiledlayout(2,1)

nexttile
plot(x,abs(dphi(x)), x,ones(size(x)),'k*')
title('derivative of dphi')

grid on

nexttile
plot(x, y, x, zeros(size(x)),'k')
grid on

            % So the plot on the buttom shows that this non-linear equation
            % has three answers around -1, 0, and 2. but the plot on the
            % top shows that the fixed point method can only give us the
            % answer which is around 0, because two other answers are in
            % the area for which the absolute value of  derivative of phi
            % function is larger than 1

x0 = 0.3;

tolerance=10.^-7;
maxIter=200;


for i = 1 : maxIter

    x1 = phi(x0);
    
    if abs( f(x1) ) < tolerance
        break
    end

    x0 = x1;
end

x1
sprintf (" the number of iterations for Newton's method = %.15g ", i)
