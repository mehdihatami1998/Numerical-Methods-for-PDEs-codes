clear all
close all
clc 
format long


ni = 0.1 %diffusity or viscosity (instead of a2)

L = 10;
N = 400;
T = 1;

alpha = 10; % for the new function

% u = 4
% 
% a1 = -u; % a1*T it says how much and in which direction we should go to take 
%         % the sampelling point. if it is negative as in this case, it will
        % go forward, and if it's positive, it will go in negative
        % direction

         
        % a1 is the opposite of velocity, and in many exercises we will be
        % asked to interpret the opposite of velocity which is u


 a2 = ni;  %the strength of this coefficient determins how fast this diffusive process
            % is tkaing place

nmodes = 100;

h = L / N;

x = [0 : h : L-h]';


% C0 is the initail condition, the new one is almost incontinous

 c0 = @(x) exp(-(x - L * 0.5) .^2 / (L * 0.1) .^2);

%c0 = @(x) -tanh(alpha*(x-L/3)) + tanh(alpha*(x - L/5));


c0hat = fft(c0(x));

% mask = ones(N,1);
% mask(nmodes+2:N-nmodes) = 0
% 
%advection equation

omega = 2 * pi / L;
%  kk = [0;[1:N/2]';[-(N/2)+1):-1]'];

kk = [0 : N/2,(-N/2+1):-1]';% this is the order that matlab makes fft
 omegak = omega * kk; % j is the imaginary unit


% dealing with zero source K
alphak = a2 *( j *  omegak).^2; %a2 and omega k are positive, j is imaginary unit
                                % it means that this square will be
                                % negative. so it will make the plot
                                % smaller


chat = ( c0hat .* exp(alphak*T));   %this is made square



% solution of time t
% c = ifft(chat.*mask);

c = ifft(chat);

plot(x, c0(x), 'b', x, c, 'r*')


% this is how the red colour spreads with time,
% the more we increase T in the first lines of code, the more it will be
% spreaded out.


% so it depends on how strong the liquid is, and how long we wait, the
% liquid will spread in the water
