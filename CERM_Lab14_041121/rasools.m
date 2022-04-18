%% Exercise 1

clc
clear all
close all
format long

A=[-30,0,-28;
   -29,-1,-29;
    0,0,-2];
g=[2,2,2]';

f=@(t,y) A*y+g;

y0=[1;2;-10];
T=4;
N=200;
h=T/N;

tn=[0:h:T]';

options=odeset('RelTol',10^-9,'MaxStep',h);
[tref,yref]=ode45(f,tn,y0,options);

plot(tref,yref(:,1),'y',tref,yref(:,2),'b',tref,yref(:,3),'r')
hold on

%bdf3 method
tic
ybdf3=zeros(size(yref));
ybdf3(1,:)=y0';
ybdf3(2,:)=yref(2,:);
ybdf3(3,:)=yref(3,:);

options=optimoptions('fsolve','Display','none','TolX',10^-8);
for k3=4:length(tn)
    uk0=ybdf3(k3-3,:)';
    uk1=ybdf3(k3-2,:)';
    uk2=ybdf3(k3-1,:)';
    
    tk3=tn(k3);
    %k es 2
    g=@(uk3) uk3-(18/11*uk2-9/11*uk1+2/11*uk0+6/11*h*f(tk3,uk3));
    
    ybdf3(k3,:)=[fsolve(g,uk2,options)]';
end

timeBDF3=toc;

plot(tn,ybdf3,'g--')
hold on

%Errors
RelErr=[];
for i=1:3
RelErr(i)=norm(yref(:,i)-ybdf3(:,i),inf)/norm(yref,inf);
end
disp('Relative error for BDF3 is max for y')
[m,pos]=max(RelErr)


%Three stage Runge Kutta (explicit method)
tic
yrk3=zeros(size(yref));
yrk3(1,:)=y0';

for k1=2:length(tn)
    uk=yrk3(k1-1,:)';
    tk=tn(k1-1);
    
    f1=f(tk,uk);
    f2=f(tk+h,uk+h*f1);
    f3=f(tk+h/2,uk+h/4*(f1+f2));
    
    yrk3(k1,:)=[uk+h/6*(f1+f2+4*f3)]';
end

tRK=toc;
plot(tn,yrk3,'c--')

%Errors
RelErr=[];
for i=1:3
RelErr(i)=norm(yref(:,i)-yrk3(:,i),inf)/norm(yref,inf);
end
disp('Relative error for Runge Kutta is max for y')
[m,pos]=max(RelErr)

disp('Duration times')
tRK
timeBDF3


%% Exercise 1b

clc
clear all
close all
format long

A=[-30,0,-28;
   -29,-1,-29;
    0,0,-2];
g=[2,2,2]';

f=@(t,y) A*y+g;

y0=[1;2;-10];
T=4;
N=20;%%%%%%%%%%%%%
h=T/N;

tn=[0:h:T]';

options=odeset('RelTol',10^-9,'MaxStep',h);
[tref,yref]=ode45(f,tn,y0,options);

plot(tref,yref(:,1),'y',tref,yref(:,2),'b',tref,yref(:,3),'r')
hold on

%bdf3 method
tic
ybdf3=zeros(size(yref));
ybdf3(1,:)=y0';
ybdf3(2,:)=yref(2,:);
ybdf3(3,:)=yref(3,:);

options=optimoptions('fsolve','Display','none','TolX',10^-8);
for k3=4:length(tn)
    uk0=ybdf3(k3-3,:)';
    uk1=ybdf3(k3-2,:)';
    uk2=ybdf3(k3-1,:)';
    
    tk3=tn(k3);
    %k es 2
    g=@(uk3) uk3-(18/11*uk2-9/11*uk1+2/11*uk0+6/11*h*f(tk3,uk3));
    
    ybdf3(k3,:)=[fsolve(g,uk2,options)]';
end

timeBDF3=toc;

plot(tn,ybdf3,'g--')
hold on

%Errors
RelErr=[];
for i=1:3
RelErr(i)=norm(yref(:,i)-ybdf3(:,i),inf)/norm(yref,inf);
end
disp('Relative error for BDF3 is max for y')
[m,pos]=max(RelErr)


%Three stage Runge Kutta (explicit method)
tic
yrk3=zeros(size(yref));
yrk3(1,:)=y0';

for k1=2:length(tn)
    uk=yrk3(k1-1,:)';
    tk=tn(k1-1);
    
    f1=f(tk,uk);
    f2=f(tk+h,uk+h*f1);
    f3=f(tk+h/2,uk+h/4*(f1+f2));
    
    yrk3(k1,:)=[uk+h/6*(f1+f2+4*f3)]';
end

tRK=toc;
plot(tn,yrk3,'c--')

%Errors
RelErr=[];
for i=1:3
RelErr(i)=norm(yref(:,i)-yrk3(:,i),inf)/norm(yref,inf);
end
disp('Relative error for Runge Kutta is max for y')
[m,pos]=max(RelErr)

disp('Duration times')
tRK
timeBDF3

%Runge Kutta fails as for N=20 results in a big h that leads to stability
%issues!

%% Exercise 2
%second order problem
%v1=v0+h*f(t0,u0,v0) QUESTION
clc
clear all
close all
format long

f=@(t,y) [y(2);-y(1)-y(2)];
y0=[1;0];
T=4;

yex=@(t) exp(-t/2).*cos(sqrt(3)/2*t)+1/sqrt(3)*exp(-t/2).*sin(sqrt(3)/2*t);

N=400;
h=T/N;
tn=[0:h:T]';

yref=yex(tn); %only y not y'
plot(tn,yref,'k')
hold on

%Heun method (average f1/f2)
tic
yHeun=zeros(length(tn),2);
yHeun(1,:)=y0';

for k1=2:length(tn)
    uk=yHeun(k1-1,:)';
    tk=tn(k1-1);
    f1=f(tk,uk);
    f2=f(tk+h,uk+h*f1);
    yHeun(k1,:)=[uk+h/2*(f1+f2)]';
end

times1(1)=toc;

plot(tn,yHeun,'g--')
hold on

%Adams Bashforth method 2step
tic
yAB2=zeros(length(tn),2);
yAB2(1,:)=y0';
yAB2(2,:)=yHeun(2,:);

for k2=3:length(tn)
    uk0=yAB2(k2-2,:)';
    tk0=tn(k2-2);
    uk1=yAB2(k2-1,:)';
    tk1=tn(k2-1);
    
    yAB2(k2,:)=[uk1+h*(3/2*f(tk1,uk1)-1/2*f(tk0,uk0))]';
end
times1(2,1)=toc;
plot(tn,yAB2,'y--')
hold on


%Theta method + BDF2
optionsfs=optimoptions('fsolve','Display','none','TolX',10^-9);

tic 
TM=0.53;
yTM=zeros(length(tn),2);
yTM(1,:)=y0';

for k1=2:length(tn)
    uk=yTM(k1-1,:)';
    tk=tn(k1-1);
    
    g=@(uk1) uk1-(uk+h*(TM*f(tk+h,uk1)+(1-TM)*f(tk,uk)));
    
    yTM(k1,:)=[fsolve(g,uk,optionsfs)]'; 
   
end

ybdf2=zeros(length(tn),2);
ybdf2(1,:)=y0';
ybdf2(2,:)=yTM(2,:);

for k2=3:length(tn)
    uk_1=ybdf2(k2-2,:)';
    uk=ybdf2(k2-1,:)';
    
    g=@(uk2) uk2-(4*uk-uk_1+2*h*f(tn(k2),uk2))/3;
    
    ybdf2(k2,:)=[fsolve(g,uk,optionsfs)]';
end

times1(3,1)=toc;

plot(tn,ybdf2,'m--')

%Second order leapfrog
tic
f=@(t,y,y1) -y-y1;  

yLeap=zeros(length(tn),2);
%yLeap1=zeros(size(tn));

yLeap(1,:)=y0'; %u0, v0
%yLeap1(1)=y0(2);
yLeap(2,1)=yLeap(1,1)+yLeap(1,2)*h; %u1=u0+v0*h
yLeap(2,2)=yLeap(1,2)+h*f(tn(1),yLeap(1,1),yLeap(1,2)); %v1=v0+h*f(t0,u0,v0) QUESTION

for k=2:length(tn)-1
    yLeap(k+1,2)=yLeap(k,2)+h*f(tn(k),yLeap(k,1),yLeap(k,2));
    
    yLeap(k+1,1)=2*yLeap(k,1)-yLeap(k-1,1)+h^2*f(tn(k),yLeap(k,1),yLeap(k,2));
end

times1(4,1)=toc

plot(tn,yLeap,'r-')

disp('Elapsed time Heun / Heun+AB_2 / Theta+BDF2 / Leapfrog')
times1

%Errors
disp('Relative errors in the whole interval Heun / Heun+AB_2 / Theta+BDF2/ Leapfrog')
RelErr=[norm(yref-yHeun(:,1),inf)/norm(yref,inf);
        norm(yref-yAB2(:,1),inf)/norm(yref,inf);
        norm(yref-ybdf2(:,1),inf)/norm(yref,inf);
        norm(yref-yLeap(:,1),inf)/norm(yref,inf)]

% Leapfrog takes the least time but gives the highest error

%% Exercise 03
%how can i include miu in the function?
clc
close all
clear all
format long

%1) Transform the equation in a 1st order system
uu=2;
f=@(t,y) [y(2);-y(1)+uu*(1-y(1)^2)*y(2)];

y0=[1;-1];
T=20;

%2)reference solution ode15s
N=500;
h=T/N;
tn=[0:h:T];

optionsod=odeset('RelTol',10^-10,'MaxStep',h);
[tref,yref]=ode15s(f,tn,y0,optionsod);
plot(tn,yref(:,1),'k',tn,yref(:,2),'k');
hold on
yref2=yref;

%3)Runge Kutta + Adams Bashforth 3 with RK

tic
yRK=zeros(length(tn),2);
yRK(1,:)=y0';

for k=1:length(tn)-1
    uk=yRK(k,:)';
    f1=f(tn(k),uk);
    f2=f(tn(k)+h,uk+h*f1);
    f3=f(tn(k)+h/2,uk+h/4*(f1+f2));
    
    yRK(k+1,:)=[uk+h/6*(f1+f2+4*f3)]';
    
end

plot(tn,yRK,'--')
hold on

yab3=zeros(length(tn),2);
yab3(1,:)=y0';
yab3(2,:)=yRK(2,:);
yab3(3,:)=yRK(3,:);

for k=3:length(tn)-1
    uk=yab3(k,:)';
    uk_1=yab3(k-1,:);
    uk_2=yab3(k-2,:);
    
    yab3(k+1,:)=uk+h/12*(23*f(tn(k),uk)-16*f(tn(k-1),uk_1)+5*f(tn(k-2),uk_2));
    
end

times1(1,1)=toc;

plot(tn,yab3,'--')
RelErr=[norm(yref(:,1)-yRK(:,1),inf)/norm(yref(:,1),inf);
        norm(yref(:,1)-yab3(:,1),inf)/norm(yref(:,1),inf)];

%Crank Nicholson + BDF2 and miu=20
uu=20;
f=@(t,y) [y(2);-y(1)+uu*(1-y(1)^2)*y(2)];
%reference solution ode15s
N=100;
h=T/N;
tn=[0:h:T];
optionsod=odeset('RelTol',10^-10,'MaxStep',h);
[tref,yref]=ode15s(f,tn,y0,optionsod);
figure()
plot(tn,yref(:,1),'k',tn,yref(:,2),'k');
hold on
yref20=yref;

tic
yCN=zeros(length(tn),2);
yCN(1,:)=y0';
optionsfs=optimoptions('fsolve','Display','none','TolX',10^-7);

for k=1:N
    uk=yCN(k,:)';
    g=@(uk1) uk1-(uk+h/2*(f(tn(k),uk)+f(tn(k+1),uk1)));
    
    yCN(k+1,:)=[fsolve(g,uk,optionsfs)]';
end

plot(tn,yCN,'-')
hold on
%stability problems with miu=20

ybdf2=zeros(length(tn),2);
ybdf2(1,:)=y0';
ybdf2(2,:)=yCN(2,:);

for k=2:N
    uk=ybdf2(k,:)';
    uk_1=ybdf2(k-1,:)';
    g=@(uk1) uk1-(4*uk-uk_1+2*h*f(tn(k+1),uk1))/3;
    ybdf2(k+1,:)=[fsolve(g,uk,optionsfs)]';
end

times1(2)=toc;
plot(tn,ybdf2,'--')
hold on
RelErr=[RelErr; norm(yref(:,1)-ybdf2(:,1),inf)/norm(yref(:,1),inf)];

%Leapfrog (Verlet)
%no need for transformation
uu=2;
f=@(t,y,y1) -y+uu*(1-y^2)*y1;
T=20;
N=500;
h=T/N;

tn=[0:h:T]';

tic
yLeap=zeros(length(tn),2);
yLeap(1,:)=y0';

%u=y v=y1
%u1=y0+y0'*h
yLeap(2,1)=yLeap(1,1)+yLeap(1,2)*h;
%v1=v0+h*f(u0,v0,t)*
yLeap(2,2)=yLeap(1,2)+h*f(tn(1),yLeap(1,1),yLeap(1,2));

for k=2:N
    %uk+1=2uk-uk-1+h^2*f(uk,vk,tk)
    yLeap(k+1,1)=2*yLeap(k,1)-yLeap(k-1,1)+h^2*f(tn(k),yLeap(k,1),yLeap(k,2));
    %vk+1=vk+h*f
    yLeap(k+1,2)=yLeap(k,2)+h*f(tn(k),yLeap(k,1),yLeap(k,2));
end

times1=[times1,toc];
figure()
%plot(tn,yref2,'k',tn,yLeap,'c--')
plot([0:T/500:T],yref2,'k')
hold on
plot(tn,yLeap,'c--')

disp('Relative error RK / AdamB3 / CN+BDF2 / Leapfrog')
RelErr=[RelErr; norm(yref2(:,1)-yLeap(:,1),inf)/norm(yref2(:,1),inf)]
    
disp('Elapsed times RK+AdamB3 / CN+BDF2 / Leapfrog(Verlet)')
times1

%Leapfrog (Verlet)
%no need for transformation
uu=20;
f=@(t,y,y1) -y+uu*(1-y^2)*y1;
T=20;
N=500;
h=T/N;

tn=[0:h:T]';

tic
yLeap=zeros(length(tn),2);
yLeap(1,:)=y0';

%u=y v=y1
%u1=y0+y0'*h
yLeap(2,1)=yLeap(1,1)+yLeap(1,2)*h;
%v1=v0+h*f(u0,v0,t)*
yLeap(2,2)=yLeap(1,2)+h*f(tn(1),yLeap(1,1),yLeap(1,2));

for k=2:N
    %uk+1=2uk-uk-1+h^2*f(uk,vk,tk)
    yLeap(k+1,1)=2*yLeap(k,1)-yLeap(k-1,1)+h^2*f(tn(k),yLeap(k,1),yLeap(k,2));
    %vk+1=vk+h*f
    yLeap(k+1,2)=yLeap(k,2)+h*f(tn(k),yLeap(k,1),yLeap(k,2));
end

times1=[times1,toc];
figure()
%plot(tn,yref2,'k',tn,yLeap,'c--')
plot([0:T/100:T],yref20,'k')
hold on
plot(tn,yLeap,'c--')

disp('Relative error RK / AdamB3 / CN+BDF2 / Leapfrog')
RelErr=[RelErr; norm(yref20(:,1)-yLeap(:,1),inf)/norm(yref20(:,1),inf)]
    
disp('Elapsed times RK+AdamB3 / CN+BDF2 / Leapfrog(Verlet)')
times1
