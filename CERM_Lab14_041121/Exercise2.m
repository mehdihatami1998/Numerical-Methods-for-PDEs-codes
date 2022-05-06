% CERM_Lab_041121
% Exercise 2
% second order equation

clc
clear all
close all
format long 

f = @(t, y) [y(2);-y(1)-y(2)];
y0 = [1;0];
T = 4;

yex=@(t) exp(-t/2) .* cos(sqrt(3) / 2 * t) + 1/sqrt(3) * ...
    exp(-t/2) .* sin(sqrt(3) / 2 * t);

N = 40;
h = T / N;

tref = [0 : h : T]';
uref = yex(tref); %This is only answer of y' doesnot include answer of y''

figure(1)
plot(tref, uref)

            % Heun Method

unHeun = zeros(N+1, 2);
tnHeun = zeros(N+1, 1);

unHeun(1,:) = y0';
tnHeun(1,1) = 0;

tic

for k = 1 : N

    tnHeun(k+1, 1) = tnHeun(k, 1) + h;

    f1 = f(tnHeun(k, 1), unHeun(k, :)');

    f2 = f(tnHeun(k+1, 1), unHeun(k, :)' + h * f1);
    
    unHeun(k+1, :)= unHeun(k, :)' + 0.5 * h * (f1 + f2);

end

timeHeun = toc
plot(tref, uref, 'b', tnHeun, unHeun, 'r*')

relErrHeun = abs(uref - unHeun(:, 1))/abs(uref);
ErrorHeun = norm(relErrHeun, inf)


            % two-step Adams Bashforth Method

                    % To provide initial conditions for this method, we
                    % can't use the reference solution because that only
                    % gives one of the columns of this function, to surve
                    % that purpose, we should use another Method which
                    % doesn't need any initial condition, Heun method for
                    % example, and use answers of that as initiating points
                    % for implementing two-step Adam Bashforth Method.


unAB2 = zeros(N+1, 2);
tnAB2 = zeros(N+1, 1);

unAB2(1,:) = y0';
tnAB2(1,1) = 0;

tnAB2(2,1) = h;
unAB2(2, :) = unHeun(2, :)';

tic
for k = 2 : N

    tnAB2(k+1, 1) = tnAB2(k, 1) + h;

    

    unAB2(k+1, :) = [unAB2(k, :)' + h * (1.5 * f(tnAB2(k, 1), unAB2(k, :)') - ...
         0.5 * f(tnAB2(k-1, 1), unAB2(k-1, :)'))]';



end
timeAdamBashforth = toc

figure(2) 

plot(tref, uref, 'b', tnAB2, unAB2, 'r--')

relErrAB2 = abs(uref - unAB2(:, 1))/abs(uref);
ErrorAB2 = norm(relErrAB2, inf)

close all

%%


            % Theta Method

theta = 0.53;

unTheta = zeros(N+1, 2);
tnTheta = zeros(N+1, 1);

unTheta(1, :) = y0';

options = optimoptions ('fsolve', 'display', 'none', 'FunctionTolerance',10^-9);

tic
for k = 1 : N

    tnTheta(k+1, 1) = tnTheta(k, 1) + h;

    g = @(y) y - unTheta(k, :)' - h * (theta * f(tnTheta(k+1, 1), y) + ...
        (1 - theta) * f(tnTheta(k, 1), unTheta(k, :)'));
    
    unTheta(k+1, :)= [fsolve(g, unTheta(k, :)', options)]';
end
timeTheta = toc

plot(tref, uref, 'b', tnTheta, unTheta, 'r*')

relErrTheta = abs(uref - unTheta(:, 1))/abs(uref);
ErrorTheta = norm(relErrTheta, inf)





            % BDF2 Method

unBDF2 = zeros(N+1, 2);
tnBDF2 = zeros(N+1, 1);

unBDF2(1, :) = y0';
unBDF2(2, :) = unTheta(2,:)';

options = optimoptions ('fsolve', 'display', 'none', 'FunctionTolerance',10^-9);

tic
for k = 2 : N

    tnBDF2(k+1, 1) = tnBDF2(k, 1) + h;

    g = @(y) y - ...
        (4/3) * unBDF2(k, :)' + (1/3) * unBDF2(k-1, :)' ...
        - (2/3) * h * f(tnBDF2(k+1, 1), y);
    
    unBDF2(k+1, :)= [fsolve(g, unBDF2(k, :)', options)]';
end

timeBDF2 = toc

figure(3)
plot(tref, uref, 'b', tnBDF2, unBDF2, 'r*')

relErrBDF2 = abs(uref - unBDF2(:, 1))/abs(uref);
ErrorBDF2 = norm(relErrBDF2, inf)




            % Leapfrog Method



unLeapfrog = zeros(N+1, 2);
tnLeapfrog = zeros(N+1, 1);

unLeapfrog(1, :) = y0'; % this is friction in physical term;
unLeapfrog(2, :) = unAB2(2, :)';
tnLeapfrog(2, :) = h;
tic
for k = 2 : N

    tnLeapfrog(k+1, 1) = tnLeapfrog(k, 1) + h;
    
    unLeapfrog(k+1, :)= [2 * unLeapfrog(k, :)' - unLeapfrog(k-1, :)' ...
        + h^2 * f(tnLeapfrog(k, 1), unLeapfrog(k, :)')]';
end

timeLeapfrog = toc

figure(4)
plot(tref, uref, 'b', tnLeapfrog, unLeapfrog, 'r*')

relErrLeapfrog = abs(uref - unLeapfrog(:, 1))/abs(uref);
ErrorLeapfrog = norm(relErrLeapfrog, inf)






            % Times

HeunMethod =          ['Time Heun Method   --> ', num2str(timeHeun), '    Error Heun Method    --> ', num2str(ErrorHeun)];
AdamBashforthMethod = ['Time AdamBashforth --> ', num2str(timeAdamBashforth), '    Error AdamBashforth  --> ', num2str(ErrorAB2 )];
ThetaMethod = ['Time Theta         --> ', num2str(timeTheta),'     Error Theta          --> ', num2str(ErrorTheta)];
BDF2Method = ['Time BDF2          --> ', num2str(timeBDF2), '     Error BDF2           --> ', num2str(ErrorBDF2)];


disp(HeunMethod)
disp(AdamBashforthMethod)
disp(ThetaMethod)
disp(BDF2Method)

 