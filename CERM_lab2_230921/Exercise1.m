%Exercise 1 after watching Lectures of the professor
close all
clear all
format long
f =  @(x) exp(-x.^ 2).* sin(2.* x + 1)
df = @(x) 2.* exp(-x.^ 2).* cos(2.*x + 1) - 2.* x.* exp(-x.^ 2).*sin(2.*x + 1)
a = -3;
b = 3;
n = 200;
h = (b-a) / n;
xmesh = [-3: h: 3]';
fxmesh = f(xmesh);
dfxmesh = zeros( size(xmesh) ); 
%makes a zeros matrix with the size of xmesh


dfxmesh(end) = ( f(b+h) - f(b) )/ h;
%because we are using forward finite different approach
%dfxmesh(end)=(f(b+h)-f(b))/h; %as for this line of code, at the point b,
%because we don't have the value for forward finite difference/ or instead of this, we can use
%backward finite difference at this pint instead.
for i = 1: n
    dfxmesh(i, 1) = ( f(xmesh(i+1, 1) ) - f(xmesh(i, 1)) ) ./ h;
end
figure(1)
plot(xmesh, df(xmesh), 'ko',xmesh, dfxmesh, 'b');
err=abs(df(xmesh)- dfxmesh);

figure(2)
semilogy(xmesh, err, 'rd-')

figure(3)
plot(xmesh, err, 'rd-')

abs_error_2 = norm (err, 2)
abs_error_inf=norm (err,inf)
rel_error_2= abs_error_2/norm(df(xmesh),2)
rel_error_inf=abs_error_inf/norm(df(xmesh),inf)

%If we want to use forward finite different method to 
%calculate a derivative, we cannot go to the last point,
%because in the calculation we have to use one step further
%and in the very last point, we can't go one step further from the last
%point.

% Examples on what are norms:
% x=[1 2 3];
% norm(x,inf) = 3;
% x=[1 2 -3];
% norm(x,inf)=3;
% norm(x,2)=3.7416...
% ukleadean: norm(x)=3.7416..
