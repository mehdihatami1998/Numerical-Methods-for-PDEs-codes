%exercise 5 CERM_Lab_200921 by me
k = [0 : 1 : 20];
h = 10.^ -k;
func = (exp(h) - 1)./ h;
abs_err = abs(1- func);
semilogy(k, abs_err, 'r*')

%because of the theory of the representaiton of floating point numbers, we
%know that numbers are rounded and they have a limited accuracy, say
%10e-16, so if our calculation and numbers go over this, the precision will
%not only increase but also decrease.
