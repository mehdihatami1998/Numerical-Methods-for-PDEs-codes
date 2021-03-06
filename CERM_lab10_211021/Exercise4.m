%%
            % Solving non linear Equation
            % Bisection Method
close all
clear all
format long

f = @(x) exp(x) - 2.* x.^2;
df = @(x) exp(x)- 4.* x;


            % Defining a mesh on the X-axis with steps of 0.01
x = [-2: 0.01: 4]';
y = f(x);

figure(1)
plot(x, y, x, zeros(size(x)),'k*')
grid on

            % The plot show there are 3 roots for this non-linear equation
            % which are around -0.5, 1.5 and 2.5, by considering x0=-0.5,
            % we'll try to find that root with bisection method in this
            % section of codes


a = -1;
b = 0;
xm = (a+b) / 2;

maxIter = 200;
tolerance = 10e-8;


while i < maxIter && abs(a-b)/ 2 > tolerance
    i=i+1;
    if f(a).*f(xm)<0
        b=xm;
    else
        a=xm;
        
    end
    xm=(a+b)/2;
end
fprintf(" Bisection Method did %.15g iterations to get x = %.15g \n" ,i , xm(end))

   
    %%
    
            % Newton Method
close all
clear all
format long
f = @(x) exp(x) - 2.* x.^2;
df = @(x) exp(x)- 4.* x;


            % Defining a mesh on the X-axis with steps of 0.01
x = [-2: 0.01: 4]';
y = f(x);

figure(1)
plot(x, y, x, zeros(size(x)),'k*')
grid on

            % The plot show there are 3 roots for this non-linear equation
            % which are around -0.5, 1.5 and 2.5, by considering x0=-0.5,
            % we'll try to find that root with Newton method in this
            % section of codes

x0 = -0.5;
maxIter = 200;
tolerance = 10e-8;

for i=1:maxIter

    delta = -f(x0) / df(x0);
    
    if abs(delta) < tolerance
        break
    end

    x0 = x0 + delta;
end

fprintf(" Newton    Method did %.15g  iterations to get x = %.15g  \n" ,i , x0)


%%
            % Secant Method
close all
clear all
format long
f = @(x) exp(x) - 2.* x.^2;
df = @(x) exp(x)- 4.* x;


            % Defining a mesh on the X-axis with steps of 0.01
x = [-2: 0.01: 4]';
y = f(x);

figure(1)
plot(x, y, x, zeros(size(x)),'k*')
grid on

            % The plot show there are 3 roots for this non-linear equation
            % which are around -0.5, 1.5 and 2.5, by considering x0=-0.5,
            % we'll try to find the root which is around -0.5, to do so, we
            % need two initial steps which are defined as below:

x1 = -0.3;
x2 = -0.8;
maxIter = 200;
tolerance= 10e-8;

for i = 1 : maxIter

    delta = -f(x2).* (x1-x2) / ( f(x1)- f(x2) );

    if abs(delta) < tolerance
        break
    end
    
    x1 = x2;
    x2 = x1+delta;
end
fprintf(" Secant    Method did %.15g  iterations to get x = %.15g  \n" ,i , x2)

%%
            %Chord Method
close all
clear all
format long

f = @(x) exp(x) - 2.* x.^2;
df = @(x) exp(x)- 4.* x;


            % Defining a mesh on the X-axis with steps of 0.01
x = [-2: 0.01: 4]';
y = f(x);

figure(1)
plot(x, y, x, zeros(size(x)),'k*')
grid on

            % The plot show there are 3 roots for this non-linear equation
            % which are around -0.5, 1.5 and 2.5, by considering x0=-0.5,
            % we'll try to find that root with Chord Method in this
            % section of codes

a= -2;
b = 4;

maxIter = 100;
tolerance = 10e-9;
x0 = 2.4;

for i = 1 : maxIter

    delta = f(x0)* (a-b) / (f(a)-f(b) );

    if abs(delta)<tolerance
        break
    else
        x0 = x0+ delta;
    end

end

fprintf(" Chord     Method did %.15g iterations to get x = %.15g  \n" ,i , x0)

%% 
            % Fixed Point Method Phi1


close all
clear all
format long
f = @(x) exp(x) - 2.* x.^2;
df = @(x) exp(x)- 4.* x;
phi1 = @(x) log( 2*x.^ 2 );
dphi1 = @(x) 2./x;


            % Defining a mesh on the X-axis with steps of 0.01
x = [-2: 0.01: 4]';
y = f(x);


figure(1)
tiledlayout(2,1)

nexttile
plot(x,abs(dphi1(x)), x,ones(size(x)),'k*')
ylim ([0 2])
title('derivative of dphi1')

grid on

nexttile
plot(x, y, x, zeros(size(x)),'k')
grid on

            % So the plot on the buttom shows that this non-linear equation
            % has three answers around -0.5, 1.5, and 2.5 but the plot on the
            % top shows that the fixed point method can only give us the
            % answer which is around 2.5, because two other answers are in
            % the area for which the absolute value of  derivative of phi
            % function is larger than 1


x0 = 2.1;
maxIter = 200;
tolerance = 10e-8;

for i = 1 : maxIter

    x1 = phi1(x0);
    if abs(x0-x1) < tolerance
        break
    end

    x0 = x1;
end

fprintf(" Fixed Point Method with Phi1 did %.15g  iterations to get x = %.15g  \n" ,i , x1)

    %%
            % Fixed Point Phi2
            close all
clear all
format long
f = @(x) exp(x) - 2.* x.^2;
df = @(x) exp(x)- 4.* x;

phi2=@(x) sqrt( 0.5* exp(x) );
dphi2=@(x) exp(x)./ ( (4.*(exp(x)./2).^ (1/2) ) );

            % Defining a mesh on the X-axis with steps of 0.01
x = [-2: 0.01: 4]';
y = f(x);


figure(1)
tiledlayout(2,1)

nexttile
plot(x,abs(dphi2(x)), x,ones(size(x)),'k*')
ylim ([0 2])
title('derivative of dphi2')

grid on

nexttile
plot(x, y, x, zeros(size(x)),'k')
grid on

            % So the plot on the buttom shows that this non-linear equation
            % has three answers around -0.5, 1.5, and 2.5 but the plot on the
            % top shows that the fixed point method can only give us the
            % answer which is around 0.5 and the one around 1.5 ,
            % because the other answer is located in 
            % the area for which the absolute value of  derivative of phi
            % function is larger than 1


x0 = 1.5;
maxIter = 200;
tolerance = 10e-8;

for i = 1 : maxIter

    x1 = phi2(x0);
    if abs(x0-x1) < tolerance
        break
    end

    x0 = x1;
end

fprintf(" Fixed Point Method with Phi2 did %.15g  iterations to get x = %.15g  \n" ,i , x1)

    %%
            % Fixed Point Phi2
            close all
clear all
format long
f = @(x) exp(x) - 2.* x.^2;
df = @(x) exp(x)- 4.* x;

phi3 = @(x) -sqrt( 0.5* exp(x) );
dphi3 = @(x) -exp(x)./ ( (4.*(exp(x)./2).^ (1/2) ) );


            % Defining a mesh on the X-axis with steps of 0.01
x = [-2: 0.01: 4]';
y = f(x);


figure(1)
tiledlayout(2,1)

nexttile
plot(x,abs(dphi3(x)), x,ones(size(x)),'k*')
ylim ([0 2])
title('derivative of dphi3')

grid on

nexttile
plot(x, y, x, zeros(size(x)),'k')
grid on

            % So the plot on the buttom shows that this non-linear equation
            % has three answers around -0.5, 1.5, and 2.5 but the plot on the
            % top shows that the fixed point method can only give us the
            % answer which is around 0.5 and the one around 1.5 ,
            % because the other answer is located in 
            % the area for which the absolute value of  derivative of phi
            % function is larger than 1


x0 = 1.5;
maxIter = 200;
tolerance = 10e-8;

for i = 1 : maxIter

    x1 = phi3(x0);
    if abs(x0-x1) < tolerance
        break
    end

    x0 = x1;
end

fprintf(" Fixed Point Method with Phi3 did %.15g  iterations to get x = %.15g  \n" ,i , x1)
