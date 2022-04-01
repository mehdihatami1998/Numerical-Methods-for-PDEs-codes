            % Solving ODEs
            % Forward Euler Method Empirical Order
clc
clear all
close all
format long

            % Defining the Problem and Exact Function
alpha = 0.001;           
f = @(t,y) (alpha * 0.5) * cos(t) * ( 1 - y.^2 );
y_exact = @(t) ( exp( alpha * sin(t) ) - 1) ./ ( exp (alpha * sin(t)) +1);


y0 = 0;
T = 6*pi;

            % Defining conditions for two solutions with h1 and h2
h1 = 0.1;
h2 = 0.05;

N1 = T/h1;
N2 = T/h2;

un1 = zeros(round(N1+1),1);
tn1 = zeros(round(N1+1),1);

un2 = zeros(size(N2+1));
tn2 = zeros(size(N2+1));

un1(1,1) = y0;
un2(1,1) = y0;

for k = 1 : N1
    un1(k+1,1) = un1(k,1) + h1 * f( tn1(k,1), un1(k,1) );
    tn1(k+1,1) = tn1(k) + h1;
end


for j = 1 : N2
    un2(j+1,1) = un2(j,1) + h2 * f( tn2(j,1), un2(j,1) );
    tn2(j+1,1) = tn2(j) + h2;
end

            % Plotting Euler Method Answers vs. Exact Answer
tmesh1 = [0 : h1 : T]';
tmesh2 = [0 : h2 : T]';

figure(1)
plot(tmesh1,y_exact(tmesh1), 'k', tn1,un1,'r*')
title("Euler Method with h = 0.1 vs. Exact answer")

figure(2)
plot(tmesh2, y_exact(tmesh2), tn2, un2, 'r*')
title("Euler Method with h=0.05 vs. Exact answer")


            % Absolute and Relatve Errors on the Interval

abs_err_EULER1 = abs( y_exact(tn1) - un1 );
abs_err_EULER2 = abs( y_exact(tn2) - un2 );


%size of rel_err_euler1 & 2 are square matrixes, not vectors

rel_err_EULER1 = abs_err_EULER1 ./ y_exact(tn1);
rel_err_EULER2 = abs_err_EULER2 

figure(3)
semilogy(tn1,rel_err_EULER1)
title("EULER Method (h=0.1) Relative Error")
figure(4)
semilogy(tn2,rel_err_EULER2)
title("EULER Method (h=0.05) Relative Error")



            % Absolute and Relative Errors at Final Time

abs_error_endpoint_EULER1 = abs( y_exact(tn1(end)) - un1(end) )
abs_error_endpoint_EULER2 = abs( y_exact(tn2(end)) - un2(end) )


rel_error_endpoint_EULER1 = abs( y_exact(tn1(end)) - un1(end)) ./ abs( y_exact( tn1(end) ) )
rel_error_endpoint_EULER2 = abs( y_exact(tn2(end)) - un2(end)) ./ abs( y_exact( tn2(end) ) )


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_EULER1 = norm( abs_err_EULER1, inf )
abs_err_Norminf_EULER2 = norm( abs_err_EULER2, inf )


rel_err_Norminf_EULER1 = norm( abs_err_EULER1, inf ) / norm( y_exact(tn1),inf )
rel_err_Norminf_EULER2 = norm( abs_err_EULER2, inf ) / norm( y_exact(tn2),inf )


            % Forward Euler Method Empirical Order Calculation 
p_emp_Euler = -log2(rel_err_Norminf_EULER2/rel_err_Norminf_EULER1)