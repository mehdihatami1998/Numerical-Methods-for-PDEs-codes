%exercise 2 
close all
clear all
format long
r = @(x) 3.* x.* log(x.^2 + 1);
dr= @(x) 3* log(x.^2 + 1) + (6.* x.^2) ./ (x.^2 + 1);
b = 1;
a = -1;
n = 50;
h = (b-a) / n;
xmesh = [a: h: b];
figure(1)
plot(xmesh, r(xmesh), 'b', xmesh, dr(xmesh), 'ro'); 
drxmesh = zeros( size(xmesh) );
for i = 2: n
    drxmesh(i) = ( r(xmesh(i+1) ) - r(xmesh(i-1)) ) ./ (2.* h);
end
figure(2)
plot(xmesh, dr(xmesh),'b', xmesh, drxmesh, 'ro');
error = (dr(xmesh)- drxmesh);
abs_err_2 = norm(error, 2);
abs_err_inf = norm(error, inf);
rel_err_2 = abs_err_2 / norm( dr(xmesh), 2 );
rel_err_inf = abs_err_inf ./ norm( dr(xmesh), inf );

    

