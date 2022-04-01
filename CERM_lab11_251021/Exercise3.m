            % Solving ODEs
            % Heun Method Empirical Order
clc
clear all
close all
format long

            % Defining the Problem and Exact Function
f = @(t,y) -y + exp(t) * sin(t);
y_exact = @(t) (51/5) * exp(-t) - 0.2 * exp(t) .* (cos(t) - 2*sin(t));

y0 = 10;
T = 5;

            % Defining conditions for two solutions with h1 and h2
h1 = 0.01;
h2 = 0.005;

N1 = T/h1;
N2 = T/h2;

un1 = zeros(round(N1+1),1);
tn1 = zeros(round(N1+1),1);

un2 = zeros(size(N2+1));
tn2 = zeros(size(N2+1));

un1(1,1) = y0;
un2(1,1) = y0;


for k = 1:N1
    tn1(k+1, 1) = tn1(k,1) + h1;
    f1 = f( tn1(k,1), un1(k,1) );
    f2 = f( tn1(k+1,1), un1(k,1)+ h1*f1);
    un1(k+1,1) = un1(k,1) + 0.5* h1 * (f1+f2);
end


for j = 1:N2
    tn2(j+1, 1) = tn2(j,1) + h2;
    f1 = f( tn2(j,1), un2(j,1) );
    f2 = f( tn2(j+1,1), un2(j,1)+ h2*f1);
    un2(j+1,1) = un2(j,1) + 0.5* h2 * (f1+f2);
end


            % Plotting Heun Method Answers vs. Exact Answer
tmesh1 = [0 : h1 : T]';
tmesh2 = [0 : h2 : T]';

figure(1)
plot(tmesh1,y_exact(tmesh1), 'k', tn1,un1,'r*')
title("Heun Method with h = 0.1 vs. Exact answer")

figure(2)
plot(tmesh2, y_exact(tmesh2), tn2, un2, 'r*')
title("Heun Method with h=0.05 vs. Exact answer")


            % Absolute and Relatve Errors on the Interval

abs_err_Heun1 = abs( y_exact(tn1) - un1 );
abs_err_Heun2 = abs( y_exact(tn2) - un2 );


    %size of rel_err_Heun1 & 2 are square matrixes, not vectors

rel_err_Heun1 = abs_err_Heun1 ./ y_exact(tn1);
rel_err_Heun2 = abs_err_Heun2 ;

figure(3)
semilogy(tn1,rel_err_Heun1)
title("Heun Method (h=0.1) Relative Error")
figure(4)
semilogy(tn2,rel_err_Heun2)
title("Heun Method (h=0.05) Relative Error")



            % Absolute and Relative Errors at Final Time

abs_error_endpoint_Heun1 = abs( y_exact(tn1(end)) - un1(end) );
abs_error_endpoint_Heun2 = abs( y_exact(tn2(end)) - un2(end) );


rel_error_endpoint_Heun1 = abs( y_exact(tn1(end)) - un1(end)) ./ abs( y_exact( tn1(end) ) );
rel_error_endpoint_Heun2 = abs( y_exact(tn2(end)) - un2(end)) ./ abs( y_exact( tn2(end) ) );


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_Heun1 = norm( abs_err_Heun1, inf );
abs_err_Norminf_Heun2 = norm( abs_err_Heun2, inf );


rel_err_Norminf_Heun1 = norm( abs_err_Heun1, inf ) / norm( y_exact(tn1),inf )
rel_err_Norminf_Heun2 = norm( abs_err_Heun2, inf ) / norm( y_exact(tn2),inf )


            % Heun Method Empirical Order Calculation 
p_emp_Heun = -log2(rel_err_Norminf_Heun2/rel_err_Norminf_Heun1)

%%
            
            % MEuler Method Empirical Order
% clc
clear all
close all
format long

            % Defining the Problem and Exact Function
f = @(t,y) -y + exp(t) * sin(t);
y_exact = @(t) (51/5) * exp(-t) - 0.2 * exp(t) .* (cos(t) - 2*sin(t));

y0 = 10;
T = 5;

            % Defining conditions for two solutions with h1 and h2
h1 = 0.01;
h2 = 0.005;

N1 = T/h1;
N2 = T/h2;

un1 = zeros(round(N1+1),1);
tn1 = zeros(round(N1+1),1);

un2 = zeros(size(N2+1));
tn2 = zeros(size(N2+1));

un1(1,1) = y0;
un2(1,1) = y0;


for k=1:N1
    tn1(k+1,1)=tn1(k,1)+h1;
    f1=f(tn1(k,1),un1(k,1));
    f2=f( tn1(k) + h1/2, un1(k,1) + 0.5 * h1 * f1);
    un1(k+1,1)=un1(k,1) + h1 * f2;
end


for k=1:N2
    tn2(k+1,1)=tn2(k,1)+h2;
    f1=f( tn2(k,1), un2(k,1));
    f2=f( tn2(k) + h2/2, un2(k,1) + 0.5 * h2 * f1);
    un2(k+1,1)=un2(k,1) + h2 * f2;
end


            % Plotting MEuler Method Answers vs. Exact Answer
tmesh1 = [0 : h1 : T]';
tmesh2 = [0 : h2 : T]';

figure(1)
plot(tmesh1,y_exact(tmesh1), 'k', tn1,un1,'r*')
title("MEuler Method with h = 0.1 vs. Exact answer")

figure(2)
plot(tmesh2, y_exact(tmesh2), tn2, un2, 'r*')
title("MEuler Method with h=0.05 vs. Exact answer")


            % Absolute and Relatve Errors on the Interval

abs_err_MEuler1 = abs( y_exact(tn1) - un1 );
abs_err_MEuler2 = abs( y_exact(tn2) - un2 );


    %size of rel_err_MEuler1 & 2 are square matrixes, not vectors

rel_err_MEuler1 = abs_err_MEuler1 ./ y_exact(tn1);
rel_err_MEuler2 = abs_err_MEuler2 ;

figure(3)
semilogy(tn1,rel_err_MEuler1)
title("MEuler Method (h=0.1) Relative Error")
figure(4)
semilogy(tn2,rel_err_MEuler2)
title("MEuler Method (h=0.05) Relative Error")



            % Absolute and Relative Errors at Final Time

abs_error_endpoint_MEuler1 = abs( y_exact(tn1(end)) - un1(end) );
abs_error_endpoint_MEuler2 = abs( y_exact(tn2(end)) - un2(end) );


rel_error_endpoint_MEuler1 = abs( y_exact(tn1(end)) - un1(end)) ./ abs( y_exact( tn1(end) ) );
rel_error_endpoint_MEuler2 = abs( y_exact(tn2(end)) - un2(end)) ./ abs( y_exact( tn2(end) ) );


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_MEuler1 = norm( abs_err_MEuler1, inf );
abs_err_Norminf_MEuler2 = norm( abs_err_MEuler2, inf );


rel_err_Norminf_MEuler1 = norm( abs_err_MEuler1, inf ) / norm( y_exact(tn1),inf )
rel_err_Norminf_MEuler2 = norm( abs_err_MEuler2, inf ) / norm( y_exact(tn2),inf )


            % MEuler Method Empirical Order Calculation 
p_emp_MEuler = -log2(rel_err_Norminf_MEuler2/rel_err_Norminf_MEuler1)


%%
            
            % leapfrog Method Empirical Order
% clc
clear all
close all
format long

            % Defining the Problem and Exact Function
f = @(t,y) -y + exp(t) * sin(t);
y_exact = @(t) (51/5) * exp(-t) - 0.2 * exp(t) .* (cos(t) - 2*sin(t));

y0 = 10;
T = 5;

            % Defining conditions for two solutions with h1 and h2
h1 = 0.01;
h2 = 0.005;

N1 = T/h1;
N2 = T/h2;

un1 = zeros(round(N1+1),1);
tn1 = zeros(round(N1+1),1);

un2 = zeros(size(N2+1));
tn2 = zeros(size(N2+1));

un1(1,1) = y0;
un2(1,1) = y0;

tn1(2,1) = h1;
tn2(2,1) = h2;

f1 = f(tn1(1,1), un1(1,1));
f2 = f(tn1(1,1) + h1/2, un1(1,1) + 0.5 * h1 * f1);
un1(2,1) = un1(1, 1) + h1 * f2;


f1 = f(tn2(1,1), un2(1,1));
f2 = f(tn2(1,1) + h2/2, un2(1,1) + 0.5 * h2 * f1);
un2(2,1) = un2(1, 1) + h2 * f2;


for k = 2 : N1
    tn1(k+1, 1) = tn1(k, 1) + h1;
    un1(k+1, 1) = un1(k-1, 1) + 2 * h1 * f(tn1(k,1),un1(k,1));
end


for j = 2 : N2
    tn2(j+1, 1) = tn2(j, 1) + h2;
    un2(j+1, 1) = un2(j-1, 1) + 2 * h2 * f(tn2(j,1),un2(j,1));
end


            % Plotting leapfrog Method Answers vs. Exact Answer
tmesh1 = [0 : h1 : T]';
tmesh2 = [0 : h2 : T]';

figure(1)
plot(tmesh1,y_exact(tmesh1), 'k', tn1,un1,'r*')
title("leapfrog Method with h = 0.1 vs. Exact answer")

figure(2)
plot(tmesh2, y_exact(tmesh2), tn2, un2, 'r*')
title("leapfrog Method with h=0.05 vs. Exact answer")


            % Absolute and Relatve Errors on the Interval

abs_err_leapfrog1 = abs( y_exact(tn1) - un1 );
abs_err_leapfrog2 = abs( y_exact(tn2) - un2 );


    %size of rel_err_leapfrog1 & 2 are square matrixes, not vectors

rel_err_leapfrog1 = abs_err_leapfrog1 ./ y_exact(tn1);
rel_err_leapfrog2 = abs_err_leapfrog2 ;

figure(3)
semilogy(tn1,rel_err_leapfrog1)
title("leapfrog Method (h=0.1) Relative Error")
figure(4)
semilogy(tn2,rel_err_leapfrog2)
title("leapfrog Method (h=0.05) Relative Error")


            % Absolute and Relative Errors at Final Time

abs_error_endpoint_leapfrog1 = abs( y_exact(tn1(end)) - un1(end) );
abs_error_endpoint_leapfrog2 = abs( y_exact(tn2(end)) - un2(end) );


rel_error_endpoint_leapfrog1 = abs( y_exact(tn1(end)) - un1(end)) ./ abs( y_exact( tn1(end) ) );
rel_error_endpoint_leapfrog2 = abs( y_exact(tn2(end)) - un2(end)) ./ abs( y_exact( tn2(end) ) );


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_leapfrog1 = norm( abs_err_leapfrog1, inf );
abs_err_Norminf_leapfrog2 = norm( abs_err_leapfrog2, inf );


rel_err_Norminf_leapfrog1 = norm( abs_err_leapfrog1, inf ) / norm( y_exact(tn1),inf )
rel_err_Norminf_leapfrog2 = norm( abs_err_leapfrog2, inf ) / norm( y_exact(tn2),inf )


            % leapfrog Method Empirical Order Calculation 
p_emp_leapfrog = -log2(rel_err_Norminf_leapfrog2/rel_err_Norminf_leapfrog1)



            % Given comparing the answers of this code, Leapfrog method is
            % the most accurate one among these three methods, then
            % Modified Euler was more accurate than Heun method, and Heun
            % method was the least accurate method here.
            %
            % Comparing the accuracy:
            % LeapFrog > Modified Euler > Heun Method
