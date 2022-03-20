
close all
clear all
format long
f = @(x) exp(-3.* x.^2);
df = @(x) - 6.* x.* exp(-3.* x.^2);

b = 2;
a = -1;

n = 150;
h = (b-a) / n;

xmesh = [a: h: b]';
figure(1)
plot(xmesh, f(xmesh));


            %forward finite difference:
dfxmesh1 = zeros( size(xmesh) );
dfxmesh1(end) = ( f(xmesh(n) + h) - f(xmesh(n)) ) ./h;

for i = 1: n - 1
    dfxmesh1(i) = (f(xmesh(i+ 1)) - f(xmesh(i)) )./ h;
end

figure(2)
plot(xmesh, df(xmesh), 'r', xmesh, dfxmesh1, 'ko')

            %backward finite difference:

dfxmesh2 = zeros(size(xmesh));
dfxmesh2(1) = ( f(xmesh(i) ) - f( xmesh(i) - h ) ) ./ h;

for i = 2: n
    dfxmesh2(i) = ( f(xmesh(i)) - f(xmesh(i-1)) ) ./ h;
end

figure(3)
plot(xmesh, df(xmesh), 'r', xmesh, dfxmesh2, 'bo')



            %centered finite difference:
dfxmesh3 = zeros(size(xmesh));

for i = 2 : n-1
    dfxmesh3(i) = ( f(xmesh(i+1)) - f(xmesh(i-1)) ) ./ (2* h);
end

figure(4)
plot(xmesh, df(xmesh), 'r', xmesh, dfxmesh3, 'ko')


            %fourth order finite difference:
dfxmesh4 = zeros(size(xmesh));
for i = 3: n-2
    dfxmesh4(i) = ( -f(xmesh(i+2)) + 8* f(xmesh(i+1)) ...
        -8* f(xmesh(i-1)) + f(xmesh(i-2)) ) ./ (12* h);
end

figure(5)
plot(xmesh, df(xmesh), 'c', xmesh, dfxmesh4, 'ro')


            %error calculations for forward finite difference

error = ( df(xmesh) - dfxmesh1 );

abs_err_2 = norm(error, 2);

abs_err_inf = norm(error, inf);

rel_err_2 = abs_err_2 / norm(df(xmesh), 2);

rel_err_inf = abs_err_inf./ norm(df(xmesh), inf);

figure(6)
semilogy(xmesh, error, 'bd-')


            %error calculation for backward finite elemes calculations:
error = (df(xmesh) - dfxmesh2);

abs_err_2 = norm(error, 2);

abs_err_inf = norm(error, inf);

rel_err_2 = abs_err_2 /norm( df(xmesh), 2 );

rel_err_inf = abs_err_inf ./ norm(df(xmesh), inf);


            %error calculations for centered finite element
error = (df(xmesh) - dfxmesh3);

abs_err_2 = norm(error, 2);

abs_err_inf = norm(error, inf);

rel_err_2 = abs_err_2 / norm( df(xmesh), 2 );

rel_err_inf = abs_err_inf ./ norm( df(xmesh), inf) ;


            %errors calculation for forth order finite difference:
error = ( df(xmesh) - dfxmesh4 );

abs_err_2 = norm(error, 2);

abs_err_inf = norm(error, inf);

rel_err_2 = abs_err_2 / norm(df(xmesh), 2);

rel_err_inf = abs_err_inf ./ norm(df(xmesh), inf);
