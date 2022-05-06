%Exercise 5
clear all
close all
format long
L=4*pi;
N=200;
h=L/N;
x=[0:h:L-h]';

y1=@(x) (2*(x.^3)-x).*cos((x.^2));
y=y1(x);
dy1=@(x) (6*(x.^2)-1).*cos(x.^2)+(2.*x).*(x-2*x.^3).*sin(x.^2);
dy=dy1(x);

yhat=fft(y);
yy1=ifft(yhat);
figure(1) %fourier series
plot(x,y,'r',x,yy1,'ko')

yhat1=fft(y)/N;
figure(2) % this is the sectrum
semilogy([1:N/2],abs(yhat(1:N/2)/N).^2,'r*')

%from this figure we can see that the spectrum is not going to Zero, so
%this is non periodic function, or even if the function itself is periodic,
%the derivative of it is non periodic

 kk=[0:N/2,-N/2+1:-1]';
omega=2.*pi/L;
yy2=ifft(omega*j*kk.*yhat);

nmodes=[10,20,40,80];
    err1=zeros(4,2);
    err=zeros(4,2);

l=1;

for i=1:length(nmodes)
    m=nmodes(i);
    mask=ones(N,1);
    mask(m+2:N-m)=0;
    
    coeff=mask.*yhat;
    FF=ifft(coeff);
    
    coeff1=(omega*j*kk).^l.*mask.*yhat;
    FF1=ifft(coeff1);
        err(i,1)=norm(y(x)-FF,inf)/norm(y(x),inf);
        err(i,2)=norm(y(x)-FF,2)/norm(y(x),2);

        err1(i,1)=norm(y1(x)-FF1,inf)/norm(y1(x),inf);
        err1(i,2)=norm(y1(x)-FF1,2)/norm(y1(x),2);
end
  pemp_2=[-log2(err(2,2)/err(1,2));-log2(err(3,2)/err(2,2));-log2(err(4,2)/err(3,2));]
    pemp_inf=[-log2(err(2,1)/err(1,1));-log2(err(3,1)/err(2,1));-log2(err(4,1)/err(3,1));]

    pemp1_2=[-log2(err1(2,2)/err1(1,2));-log2(err1(3,2)/err1(2,2));-log2(err1(4,2)/err1(3,2));]
    pemp1_inf=[-log2(err1(2,1)/err1(1,1));-log2(err1(3,1)/err1(2,1));-log2(err1(4,1)/err1(3,1));]
