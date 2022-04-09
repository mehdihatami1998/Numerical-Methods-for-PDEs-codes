            % Solving ODEs
            % Using ODE45 Solver
% clc
clear all
close all
format long

            % Defining the Problem and Exact Function
f = @(t,y) -y + exp(t) * sin(t);
yex = @(t) (51/5) * exp(-t) - 0.2 * exp(t) .* (cos(t) - 2*sin(t));

y0 = 10;
T = 5;


options_ode45 = odeset('RelTol', 10e-7, 'MaxStep',0.01);
[tn, un] = ode45 (f, [0: 0.01: 5], 10, options_ode45);

tmesh = [0: 0.01: T]';
            % Plotting ODE45 Method approximation vs the Exact answer
figure(5)
tmesh = [0: 0.01: T]';
plot(tmesh,yex(tmesh),'r',tn,un,'ko')
title("ODE45 Method vs. Exact function")

            % Absolute and Relatve Errors on the Interval

abs_err_ODE45 = abs( yex(tn) - un );
rel_err_ODE45 = abs_err_ODE45 ./ yex(tn);

figure(6)
semilogy(tn,rel_err_ODE45)
title("ODE45 Method , Relative Error")

            % Absolute and Relative Errors at Final Time

abs_error_endpoint_ODE45 = abs( yex(tn(end)) - un(end) )
rel_error_endpoint_ODE45 = abs( yex(tn(end)) - un(end)) / abs( yex( tn(end) ) )


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_ODE45 = norm( abs_err_ODE45, inf )
rel_err_Norminf_ODE45 = norm( abs_err_ODE45, inf ) / norm( yex(tn))

%%
            % Same answer without giving Maximum time step to the 
            % ODE45 solver

            % Solving ODEs
            % Using ODE45 Solver
% clc
clear all
close all
format long

            % Defining the Problem and Exact Function
f = @(t,y) -y + exp(t) * sin(t);
yex = @(t) (51/5) * exp(-t) - 0.2 * exp(t) .* (cos(t) - 2*sin(t));

y0 = 10;
T = 5;


options_ode45 = odeset('RelTol', 10e-7); %from here I have removed 
[tn, un] = ode45 (f, [0: 0.01: 5], 10, options_ode45);

tmesh = [0: 0.01: T]';
            % Plotting ODE45 Method approximation vs the Exact answer
figure(5)
tmesh = [0: 0.01: T]';
plot(tmesh,yex(tmesh),'r',tn,un,'ko')
title("ODE45 Method vs. Exact function")

            % Absolute and Relatve Errors on the Interval

abs_err_ODE45 = abs( yex(tn) - un );
rel_err_ODE45 = abs_err_ODE45 ./ yex(tn);

figure(6)
semilogy(tn,rel_err_ODE45)
title("ODE45 Method , Relative Error")

            % Absolute and Relative Errors at Final Time

abs_error_endpoint_ODE45 = abs( yex(tn(end)) - un(end) )
rel_error_endpoint_ODE45 = abs( yex(tn(end)) - un(end)) / abs( yex( tn(end) ) )


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_ODE45 = norm( abs_err_ODE45, inf )
rel_err_Norminf_ODE45 = norm( abs_err_ODE45, inf ) / norm( yex(tn))
