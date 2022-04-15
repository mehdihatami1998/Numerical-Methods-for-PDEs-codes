            % Solving non-linear equations
            % Fixed-Point Method
clear all
close all
format long

f = @(x) (2* x.^ 2- 3* x- 2) / (x-1);
df = @(x) ( - 2* x.^2 + 3*x + 2) ./ (x - 1).^2 + (4.*x - 3)./(x - 1);
phi1 = @(x) x - 2 + ( x/(x-1) );
dphi1 = @(x) 1./ (x - 1) - x./ (x - 1).^2 + 1;

            % Defining a mesh on the X-axis with steps of 0.01

x=[-3: 0.01: 3]';
y=f(x);

tiledlayout(2,1)
nexttile

plot(x, abs(dphi1(x)), x, ones(size(x)),'k*')
ylim ([-1 2])
grid on
nexttile

plot(x, y, x, zeros(size(x)),'k')
grid on
            % apprximate answers based on the plot are: -0.5 , 2

tolerance=10e-8;
maxIter=200;
x0=1.2;
for i=1:maxIter
    i
    x1=phi1(x0);
    if abs(x0-x1) / x0 < tolerance
        break
    end
    x0=x1;
end
fprintf(" Fixed Point Method with Phi 1 did %.15g  iterations to get x = %.15g  \n" ,i , x0)

%%
            % Fixed-point Method with Phi 2
            % This fixed-point formula doesn't seem to be working

phi2 = @(x) (3 .* x.^2 - 5.*x) / (x-1);
dphi2 = @(x) ((6*x-5).*(x-1)-((3 * x.^2 - 5*x)))/((x-1).^2);
            % Defining a mesh on the X-axis with steps of 0.01

            
tiledlayout(2,1)
nexttile

plot(x, abs(dphi2(x)), x, ones(size(x)),'k*')
ylim ([-1 2])
grid on
nexttile

plot(x, y, x, zeros(size(x)),'k')
grid on
            % apprximate answers based on the plot are: -0.5 , 2

tolerance=10e-8;
maxIter=200;
x0=-0.8;
for i=1:maxIter
    i
    x1=phi2(x0);
    if abs(x0-x1)<tolerance
        break
    end
    x0=x1;
end
fprintf(" Fixed Point Method with Phi 2 did %.15g  iterations to get x = %.15g  \n" ,i , x0)


