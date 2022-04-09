            % Solving ODEs
            % Using ODE45 Solver
clc
clear all
close all
format long

            % Defining the Problem and Exact Funtion
f = @(t,y) -2 * t * y;
yex = @(t) exp(-(t.^2)) / 2;

y0 = 1/2;
T = 2;


options_ode45 = odeset('RelTol', 10e-10, 'MaxStep', 0.01);
[tn_ode45, un_ode45] = ode45 (f, [0: 0.01: 2], 0.5, options_ode45);


            % Plotting ODE45 Method approximation vs the Exact answer
tmesh = [0: 0.01: T]';
figure(1)
tmesh = [0: 0.01: T]';
plot(tmesh,yex(tmesh),'r',tn_ode45,un_ode45,'ko')
title("ODE45 Method vs. Exact function")

            % Absolute and Relatve Errors on the Interval

abs_err_ODE45 = abs( yex(tn_ode45) - un_ode45 );
rel_err_ODE45 = abs_err_ODE45 ./ yex(tn_ode45);

figure(2)
semilogy(tn_ode45,rel_err_ODE45)
title("ODE45 Method , Relative Error")

            % Absolute and Relative Errors at Final Time

abs_error_endpoint_ODE45 = abs( yex(tn_ode45(end)) - un_ode45(end) );
rel_error_endpoint_ODE45 = abs( yex(tn_ode45(end)) - un_ode45(end)) / abs( yex( tn_ode45(end) ) );


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_ODE45 = norm( abs_err_ODE45, inf );
rel_err_Norminf_ODE45 = norm( abs_err_ODE45, inf ) / norm( yex(tn_ode45))

%%
            % Solving ODEs
            % Using ode23 Solver
% clc
clear all
close all
format long

            % Defining the Problem and Exact Function
f = @(t,y) -2 * t * y;
yex = @(t) exp(-(t.^2)) / 2;

y0 = 1/2;
T = 2;


options_ode23 = odeset('RelTol', 10e-10,'MaxStep',0.1);
[tn_ode23, un_ode23] = ode23 (f, [0: 0.1: 2], 0.5, options_ode23);


            % Plotting ode23 Method approximation vs the Exact answer
tmesh = [0: 0.01: T]';
figure(1)
tmesh = [0: 0.01: T]';
plot(tmesh,yex(tmesh),'r',tn_ode23,un_ode23,'ko')
title("ode23 Method vs. Exact function")

            % Absolute and Relatve Errors on the Interval

abs_err_ode23 = abs( yex(tn_ode23) - un_ode23 );
rel_err_ode23 = abs_err_ode23 ./ yex(tn_ode23);

figure(2)
semilogy(tn_ode23,rel_err_ode23)
title("ode23 Method , Relative Error")

            % Absolute and Relative Errors at Final Time

abs_error_endpoint_ode23 = abs( yex(tn_ode23(end)) - un_ode23(end) );
rel_error_endpoint_ode23 = abs( yex(tn_ode23(end)) - un_ode23(end)) / abs( yex( tn_ode23(end) ) );


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_ode23 = norm( abs_err_ode23, inf );
rel_err_Norminf_ode23 = norm( abs_err_ode23, inf ) / norm( yex(tn_ode23))

%%

            % Solving ODEs
            % Using ode113 Solver
% clc
clear all
close all
format long

            % Defining the Problem and Exact Fun_ode113ction
f = @(t,y) -2 * t * y;
yex = @(t) exp(-(t.^2)) / 2;

y0 = 1/2;
T = 2;


options_ode113 = odeset('RelTol', 10e-10);
[tn_ode113, un_ode113] = ode113 (f, [0: 0.01: 2], 0.5, options_ode113);


            % Plotting ode113 Method approximation vs the Exact answer
tmesh = [0: 0.01: T]';
figure(1)
tmesh = [0: 0.01: T]';
plot(tmesh,yex(tmesh),'r',tn_ode113,un_ode113,'ko')
title("ode113 Method vs. Exact function")

            % Absolute and Relatve Errors on the Interval

abs_err_ode113 = abs( yex(tn_ode113) - un_ode113 );
rel_err_ode113 = abs_err_ode113 ./ yex(tn_ode113);

figure(2)
semilogy(tn_ode113,rel_err_ode113)
title("ode113 Method , Relative Error")

            % Absolute and Relative Errors at Final Time

abs_error_endpoint_ode113 = abs( yex(tn_ode113(end)) - un_ode113(end) );
rel_error_endpoint_ode113 = abs( yex(tn_ode113(end)) - un_ode113(end)) / abs( yex( tn_ode113(end) ) );


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_ode113 = norm( abs_err_ode113, inf );
rel_err_Norminf_ode113 = norm( abs_err_ode113, inf ) / norm( yex(tn_ode113))

%%
            % Solving ODEs
            % Using ode15s Solver
% clc
clear all
close all
format long

            % Defining the Problem and Exact Fun_ode15sction
f = @(t,y) -2 * t * y;
yex = @(t) exp(-(t.^2)) / 2;

y0 = 1/2;
T = 2;


options_ode15s = odeset('RelTol', 10e-10);
[tn_ode15s, un_ode15s] = ode15s (f, [0: 0.01: 2], 0.5, options_ode15s);


            % Plotting ode15s Method approximation vs the Exact answer
tmesh = [0: 0.01: T]';
figure(1)
tmesh = [0: 0.01: T]';
plot(tmesh,yex(tmesh),'r',tn_ode15s,un_ode15s,'ko')
title("ode15s Method vs. Exact function")

            % Absolute and Relatve Errors on the Interval

abs_err_ode15s = abs( yex(tn_ode15s) - un_ode15s );
rel_err_ode15s = abs_err_ode15s ./ yex(tn_ode15s);

figure(2)
semilogy(tn_ode15s,rel_err_ode15s)
title("ode15s Method , Relative Error")

            % Absolute and Relative Errors at Final Time

abs_error_endpoint_ode15s = abs( yex(tn_ode15s(end)) - un_ode15s(end) );
rel_error_endpoint_ode15s = abs( yex(tn_ode15s(end)) - un_ode15s(end)) / abs( yex( tn_ode15s(end) ) );


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_ode15s = norm( abs_err_ode15s, inf );
rel_err_Norminf_ode15s = norm( abs_err_ode15s, inf ) / norm( yex(tn_ode15s))

%%
            % Solving ODEs
            % Using ode23tb Solver
% clc
clear all
close all
format long

            % Defining the Problem and Exact Fun_ode23tbction
f = @(t,y) -2 * t * y;
yex = @(t) exp(-(t.^2)) / 2;

y0 = 1/2;
T = 2;


options_ode23tb = odeset('RelTol', 10e-10);
[tn_ode23tb, un_ode23tb] = ode23tb (f, [0: 0.01: 2], 0.5, options_ode23tb);


            % Plotting ode23tb Method approximation vs the Exact answer
tmesh = [0: 0.01: T]';
figure(1)
tmesh = [0: 0.01: T]';
plot(tmesh,yex(tmesh),'r',tn_ode23tb,un_ode23tb,'ko')
title("ode23tb Method vs. Exact function")

            % Absolute and Relatve Errors on the Interval

abs_err_ode23tb = abs( yex(tn_ode23tb) - un_ode23tb );
rel_err_ode23tb = abs_err_ode23tb ./ yex(tn_ode23tb);

figure(2)
semilogy(tn_ode23tb,rel_err_ode23tb)
title("ode23tb Method , Relative Error")

            % Absolute and Relative Errors at Final Time

abs_error_endpoint_ode23tb = abs( yex(tn_ode23tb(end)) - un_ode23tb(end) );
rel_error_endpoint_ode23tb = abs( yex(tn_ode23tb(end)) - un_ode23tb(end)) / abs( yex( tn_ode23tb(end) ) );


            % Absolute and Relative Maximum Norm on the Interval

abs_err_Norminf_ode23tb = norm( abs_err_ode23tb, inf );
rel_err_Norminf_ode23tb = norm( abs_err_ode23tb, inf ) / norm( yex(tn_ode23tb))
