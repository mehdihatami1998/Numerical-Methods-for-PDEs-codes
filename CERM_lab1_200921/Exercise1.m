%  Exercise_1 Lab 20 sept 2021 thursday
f= @(a,x) a.*exp(x)
f_approx= @(a, x) a.* (1+ x+ x.^2 / 2)
F_approx_second= @(a, x) a.* (1+ x+ x^2/2+ x^3/6+ x^4/24)
%
%1-1
abs_err= f(2, 0.2) - f_approx(2, 0.2)
rel_err=(f(2, 0.2)- f_approx(2, 0.2))./ f(2, 0.2)
%
%1-2
abs_err=f(10^5, 0.2) - f_approx(10^5, 0.2)
rel_err=(f(10^5, 0.2)-f_approx(10^5, 0.2)) ./ f(10^5, 0.2)
%
%1-3
abs_err=f(10^3, 4) - f_approx(10^3, 4)
rel_err=(f(10^3,4) - f_approx(10^3, 4))./ f(10^3, 4)
%1-4
%1-4-1
abs_err=f(2, 0.2) - F_approx_second(2, 0.2)
rel_err=(f(2, 0.2)-F_approx_second(2, 0.2))./ f(2, 0.2)
%1-4-2
abs_err=f(10^5, 0.2) - F_approx_second(10^5, 0.2)
rel_err=(f(10^5, 0.2)- F_approx_second(10^5, 0.2))./ f(10^5, 0.2)
%1-4-3
abs_err=f(10^3, 4) - F_approx_second(10^3, 4)
rel_err=(f(10^3, 4)- F_approx_second(10^3, 4))./f(10^3, 4)