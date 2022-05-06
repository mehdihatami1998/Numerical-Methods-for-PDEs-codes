% CERM_Lab_151121
% Exercise 2

clc
clear all
close all
format long



%%%%%%%%%%%%%%%% Comments %%%%%%%%%%%%%%%%
% If we use fourier series to approximate a regular function, the error
% will go to zero extremely fast, so fast that there won't be a constant
% number for the value of empirical estimate of the convergence order 



N1 = 100;
N2 = 50;
L = 10;
h1 = L / N1;
h2 = L / N2;
m1 = 24; 

nmode1 = m1;

x1 = [0 : h1 : L-h1]';
x2 = [0 : h2 : L-h2]';



            % Keep one of these pair of functions uncommented to calculate
            % on, and keep others commented to be neglected by the program.

% Pair 1
% y1 = exp(-abs((x1 - L/2) / (L/20)) );
% y1_2=exp(-abs((x2 - L/2) / (L/20)) );


% Pair 2
y1 = exp( -((x1 - L/2) / (L/10)).^2);
y1_2 = exp( -((x2 - L/2) / (L/10)).^2);


% Pair 3
% y1 = exp( -abs((x1 - L/2) / (L/20)) );
% y1_2 = exp( -abs((x2 - L/2) / (L/20)) );


                    % Comment
                    % here that we are adding noise to the function, it is
                    % not derivable anymore 



% Pair 4
% y1 = exp( -((x1 - L/2) / (L/10)).^2) + 0.001 * randn(N1, 1);
% y1_2 = exp( -((x2 - L/2) / (L/10)).^2) + 0.001 * randn(N2, 1);


% % Pair 5
% y1 = exp( -((x1 - L/2) / (L/10)) .^2) + 0.1 * randn(N1 ,1);
% y1_2 = exp( -((x2 - L/2) / (L/10)) .^2) + 0.1 * randn(N2 ,1);

yhat1=fft(y1);
yhat1_2=fft(y1_2);

figure(1)
plot([1:N1/2],abs(yhat1(1:N1/2)).^2,'bo')

mask1=ones(N1,1);
mask1(nmode1+2:N1-(nmode1))=0;
mask_2=ones(N2,1);
mask_2(nmode1+2:N2-(nmode1))=0;



yy1=ifft(yhat1.*mask1);
yy1_2=ifft(yhat1_2.*mask_2);


norm_two_1=norm(y1-yy1,2);
norm_inf_1=norm(y1-yy1,inf)


norm_two_1_2=norm(y1_2-yy1_2,2);
norm_inf_1_2=norm(y1_2-yy1_2,inf)


p_emp_y1=-log2(norm_inf_1_2/norm_inf_1)


