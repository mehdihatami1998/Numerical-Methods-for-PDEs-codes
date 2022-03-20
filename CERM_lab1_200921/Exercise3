%using a for loop
n=100;
X=zeros(1,n)';
E=zeros(1,n)';
a=2;
X(1)=0.5;
tolerance=1e-15;
for i=1:100
    X(i+1)= 0.5*(X(i) + a/X(i));
    E(i+1) =abs(X(i+1) - X(i));
    if E(i+1)<tolerance
        disp('done');
        break
    end
end
X(i:end)=[];
X
%%
%using a while loop
%
clear all
close all
clc
n=100;
X=zeros(1,n)';
E=zeros(1,n)';
a=2;
X(1)=0.5;
tolerance=1e-20;
tol=1;
i=1;
while i<100
    X(i+1)= 0.5*(X(i) + a/X(i));
    E(i+1) =abs(X(i+1) - X(i));
    if E(i+1)<tolerance
        disp('done');
        break
    end
    i=i+1
end
X(i:end)=[];
X
