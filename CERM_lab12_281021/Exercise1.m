% CERM_Lab_281021
% Exercise 1
            
            % Solving ODEs
            % Three-step Runge-Kutta (RK) Mehtod

clc
clear all
close all
format long

f = @(t,y) -y + exp(t) * sin(t);
yex = @(t) (51/5) * exp(-t) - 0.2 * exp(t) .* (cos(t) - 2 * sin(t));

y0 = 10;
T = 5;

            % Defining conditions for two solutions with N1 and N2

N1 = 5000; % This is considered RK1
N2 = 2500; % This is considered RK2

h1 = T/N1;
h2 = T/N2;

un1 = zeros(round(N1+1),1);
tn1 = zeros(round(N1+1),1);

un2 = zeros(round(N2+1),1);
tn2 = zeros(round(N2+1),1);

un1(1,1) = y0;
un2(1,1) = y0;

tn1(2, 1) = h1;
tn2(2, 1) = h2;

un1(2,1) = yex(tn1(2, 1));
un2(2,1) = yex(tn2(2, 1));

tic
for k = 2 : N1
    tn1(k+1, 1) = tn1(k, 1) + h1;
    f1 = f(tn1(k, 1), un1(k, 1));
    f2 = f(tn1(k, 1) + h1, un1(k, 1) + h1 * f1);
    f3 = f(tn1(k, 1) + h1 / 2 , un1(k, 1) + 0.25 * h1 * (f1 + f2));
    un1(k+1, 1) = un1(k, 1) + h1 * (f1 + f2 + 4 * f3) / 6;
end


for j = 2 : N2
    tn2(j+1, 1) = tn2(j, 1) + h2;
    f1 = f(tn2(j, 1), un2(j, 1));
    f2 = f(tn2(j, 1) + h2, un2(j, 1) + h2 * f1);
    f3 = f(tn2(j, 1) + h2/2 ,un2(j, 1) + 0.25 * h2 * (f1 + f2));
    un2(j+1, 1) = un2(j, 1) + h2 * (f1 + f2 + 4 * f3) / 6;
end
toc
            % Plotting RK Method Answers vs. Exact Answer
tmesh1 = [0 : h1 : T]';
tmesh2 = [0 : h2 : T]';

figure(1)
plot(tmesh1,yex(tmesh1), 'k', tn1,un1,'r*')
title("RK Method with N=5000 vs. Exact answer")

figure(2)
plot(tmesh2, yex(tmesh2), tn2, un2, 'r*')
title("RK Method with N=2500 vs. Exact answer")


            % Absolute and Relatve Errors on the Interval

abs_err_RK1 = abs( yex(tn1) - un1 );
abs_err_RK2 = abs( yex(tn2) - un2 );



rel_err_RK1 = abs_err_RK1 ./ yex(tn1);
rel_err_RK2 = abs_err_RK2 ./ yex(tn2);


figure(3)
semilogy(tn1,rel_err_RK1)
title("RK Method (N=5000) Relative Error")
figure(4)
semilogy(tn2,rel_err_RK2)
title("RK Method (N=2500) Relative Error")



            % Absolute and Relative Errors at Final Time

abs_error_endpoint_RK1 = abs( yex(tn1(end)) - un1(end) );
abs_error_endpoint_RK2 = abs( yex(tn2(end)) - un2(end) );


rel_error_endpoint_RK1 = abs( yex(tn1(end)) - un1(end)) ./ abs( yex( tn1(end) ) );
rel_error_endpoint_RK2 = abs( yex(tn2(end)) - un2(end)) ./ abs( yex( tn2(end) ) );


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_RK1 = norm( abs_err_RK1, inf );
abs_err_Norminf_RK2 = norm( abs_err_RK2, inf );


rel_err_Norminf_RK1 = norm( abs_err_RK1, inf ) / norm( yex(tn1),inf );
rel_err_Norminf_RK2 = norm( abs_err_RK2, inf ) / norm( yex(tn2),inf );


            % Forward RK Method Empirical Order Calculation 
p_emp_RK3 = -log2(rel_err_Norminf_RK1/rel_err_Norminf_RK2)

%%
            
            % Solving ODEs
            % Three-step Adams-Bashforth (AB) Mehtod


clear all
close all
format long

f = @(t,y) -y + exp(t) * sin(t);
yex = @(t) (51/5) * exp(-t) - 0.2 * exp(t) .* (cos(t) - 2 * sin(t));

y0 = 10;
T = 5;

            % Defining conditions for two solutions with N1 and N2

N1 = 5000; % This is considered AB1
N2 = 2500; % This is considered AB2

h1 = T/N1;
h2 = T/N2;

un1 = zeros(round(N1+1),1);
tn1 = zeros(round(N1+1),1);

un2 = zeros(round(N2+1),1);
tn2 = zeros(round(N2+1),1);

un1(1,1) = y0;
un2(1,1) = y0;

tn1(2, 1) = h1;
tn2(2, 1) = h2;

un1(2,1) = yex(tn1(2, 1));
un2(2,1) = yex(tn2(2, 1));

tn1(3, 1) = 2 * h1;
tn2(3, 1) = 2 * h2;

un1(3,1) = yex(tn1(3, 1));
un2(3,1) = yex(tn2(3, 1));


tic
for k = 3 : N1
    tn1(k+1, 1) = tn1(k, 1) + h1;
    un1(k+1) = un1(k, 1) + h1 * ...
        (23/12 * f(tn1(k, 1),un1(k, 1)) ...
        - 16/12 * f(tn1(k - 1), un1(k - 1))  ...
        + 5/12 * f(tn1(k - 2), un1(k - 2)));
end


for j = 3 : N2
    tn2(j+1, 1) = tn2(j, 1) + h2;
    un2(j+1) = un2(j, 1) + h2 *[ ...
        23/12 * f(tn2(j, 1),un2(j, 1)) ...
        - 16/12 * f(tn2(j - 1), un2(j - 1))  ...
        + 5/12 *f(tn2(j - 2), un2(j - 2)) ] ;
    
end
toc
            % Plotting AB Method Answers vs. Exact Answer
tmesh1 = [0 : h1 : T]';
tmesh2 = [0 : h2 : T]';

figure(1)
plot(tmesh1,yex(tmesh1), 'k', tn1,un1,'r*')
title("AB Method with N=5000 vs. Exact answer")

figure(2)
plot(tmesh2, yex(tmesh2), tn2, un2, 'r*')
title("AB Method with N=2500 vs. Exact answer")


            % Absolute and Relatve Errors on the Interval

abs_err_AB1 = abs( yex(tn1) - un1 );
abs_err_AB2 = abs( yex(tn2) - un2 );



rel_err_AB1 = abs_err_AB1 ./ yex(tn1);
rel_err_AB2 = abs_err_AB2 ./ yex(tn2);


figure(3)
semilogy(tn1,rel_err_AB1)
title("AB Method (N=5000) Relative Error")
figure(4)
semilogy(tn2,rel_err_AB2)
title("AB Method (N=2500) Relative Error")



            % Absolute and Relative Errors at Final Time

abs_error_endpoint_AB1 = abs( yex(tn1(end)) - un1(end) );
abs_error_endpoint_AB2 = abs( yex(tn2(end)) - un2(end) );


rel_error_endpoint_AB1 = abs( yex(tn1(end)) - un1(end)) ./ abs( yex( tn1(end) ) );
rel_error_endpoint_AB2 = abs( yex(tn2(end)) - un2(end)) ./ abs( yex( tn2(end) ) );


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_AB1 = norm( abs_err_AB1, inf );
abs_err_Norminf_AB2 = norm( abs_err_AB2, inf );


rel_err_Norminf_AB1 = norm( abs_err_AB1, inf ) / norm( yex(tn1),inf );
rel_err_Norminf_AB2 = norm( abs_err_AB2, inf ) / norm( yex(tn2),inf );


            % Forward AB Method Empirical Order Calculation 
p_emp_AB3 = -log2(rel_err_Norminf_AB1/rel_err_Norminf_AB2)
close all