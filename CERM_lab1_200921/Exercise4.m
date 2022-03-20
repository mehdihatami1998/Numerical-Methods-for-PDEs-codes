
n = 100;
x = zeros(n,1);
e = zeros(n,1);
a = 16;
x(1) = 0.4;
for i = 1:n
    x(i+1) = 0.5.* x(i) .* (3-a.* x(i).^2);
    if abs( x(i+1) - x(i) )/ abs( x(i+1) ) < eps* 10^3
        break
    end
end
x(i:end) = [];
x