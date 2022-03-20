 %Exercise_2 Lab 20 sept 2021 thursday
fun_abs=@(k) (1+ 10.^ (-k)- 1);
fun_approx=@(k) (abs(((1+ 10^ -k)-1)- 10^-k)) / 10^-k;
%
%
%version 1 with using loops
err_vector=zeros(1, 20)';
key=zeros(1, 20)';
for k = 1 : 20
    rel_err= abs((fun_abs(k)- fun_approx(k)))/ abs(fun_abs(k));
   err_vector(k)=rel_err
    key(k)=k
end

figure(2) 
semilogy(key,err_vector, 'r*')

%new code without loops only using vectors

k=[0:1:20]';
func1= (1+10.^(-k))-1;
approx2= abs(((1+10.^(-k)-1)-10.^(-k)))./(10.^(-k));
relative_error1= abs(func1-approx2)./abs(func1);
figure(1)
semilogy(k,relative_error1,'ro');



% part 2 of exercise 2


k=[0:1:20]';
func1= (1+10.^(-k))-1;
approx2= abs(1+10.^(-k)-10.^(-k)-1)./(10.^(-k));
relative_error2= abs(func1-approx2)./abs(func1);
figure(2)= semilogy(k,relative_error2,'ro');


