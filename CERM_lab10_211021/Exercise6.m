            % Non-linear equation
            % Fixed-point method
clear all
close all
clear all

f = @(x) sin(x)+ 0.1- 2.* x+ x.^3;
df = @(x) cos(x) - 2 + 3.*x.^2;
phi = @(x) 0.5.* ( sin(x) + x.^3 + 0.1 );
dphi = @(x) cos(x)./2 + (3.*x.^2)/2;

x = [0: 0.01: 2]';
y = f(x);

tiledlayout(2,1)
nexttile

plot(x, abs(dphi(x)), x, ones(size(x)),'k*')
ylim ([-1 2])
grid on
nexttile

plot(x, y, x, zeros(size(x)),'k')
grid on

            % There are two answers for this eqations, but we can only find
            % the one which is around 0, based on the theory, fixed-point
            % method can't give us the other answer.

x0 = 0.5;
tolerance = 10e-9;
maxIter=200;

for i = 1 : maxIter
    x1=phi(x0);

    if abs(x1-x0) < tolerance
        break
    end

    x0 = x1;
end
fprintf(" Fixed Point Method with Phi 1 did %.15g  iterations to get x = %.15g  \n" ,i , x1)

%%
            % Newton Method
clear all
close all
format long

f = @(x) sin(x) + 0.1- 2.* x+ x.^3;
df = @(x) cos(x) - 2 + 3.*x.^2;


x = [0: 0.01: 2]';
y = f(x);

figure(1)
plot(x, y, x, zeros(size(x)), 'k')
grid on

            % Given the printed graph, there are 2 answers for this
            % equation, one around 0.1, and the other around 1, now let's
            % considre x0 = 0 to find the first answer

x0 = 0;
tolerance = 10e-9;
maxIter = 200;

for i = 1 : maxIter

    delta= -f(x0) / df(x0);

    if abs(delta) < tolerance
        break
    end
    x0 = x0 + delta;

end

fprintf(" Newton Method did %.15g  iterations to get x = %.15g  \n" ,i , x0)


            % Using fsolve function to solve the equation


options=optimset('MaxIter',1000, 'TolFun',10e-11);

fsolve(f, 0, options)