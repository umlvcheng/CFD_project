
%% linerized shallow water equation

clear;
close all;
clc;
% standing wave case

% define variable
depth = 10;
a = 1280; % longth of the bay
b = 640;  % width of the bay
alpha = 0.1;
m = 2;
n = 1;
g = 10; % m/s^2;


% define the domain

xmin = -a/2; xmax = a/2;
ymin = -b/2; ymax = b/2;

% define the BC of the linerized SWE
periodic = 1;

% define the compute grid

dx = 8; dy = 8; % grid spacing
M = a/dx; % cell number at x direction
N = b/dy; % cell number at y direction

xe = xmin:dx:xmax; % cell edge at x direction
ye = ymin:dy:ymax; % cell edge at y direction

xc = (xmin+dx/2):dx:(xmax - dx/2); % eta location x
yc = ymin+dy/2:dy:ymax - dy/2; % eta locatoin y

xev = xc; % v location x
yev = ye; % v locaiton y

am = [0.5  0.5];      % averaging coefficients
ad = [1/dx -1/dx];      % gradient coefficients

% Define integration time and number of time steps
H = 10;
depth = H*ones([M N]);
Courant  = 0.5;              % Set courant number
c  = sqrt( g*max( depth(:)));              % maximum gravity wave speed
dt = Courant*dx/c;           % Courant number is fixed
ntimesteps=2*(xmax-xmin)/c/dt; % enough steps for a SINGLE circuit

% initial condition

ll = 25; % const

u = zeros([M+1 N]);
v = zeros([M N+1]);

for i = 1:length(xc)
    for j = 1:length(yc)
        rr = (xc(i) - 0)^2 + (yc(j) - 0)^2;
        eta(i,j) = alpha*exp((-1)*rr/ll^2);
    end
end

fcoriolis = zeros([M+1 N+1]);

%%  shallow water equation

for n = 1:ntimesteps;
    
   [u,v,eta] = swerk3(u,v,eta, dx,dy,dt, g, depth, fcoriolis);
    
%     if n * dt == 32 ;
%         sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
%         pause();
%     elseif n * dt == 64 ;
%         sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
%         pause();
%     elseif n * dt == 128 ;
%         sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
%         pause();
%     elseif n * dt == 256 ;
%         sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
%         pause();
%         
%     end
    % monitor the total volume of displace water
    Vol(n) = dx*dy*sum(eta(:));
    
    % monitor the energy defined on the p=points
    En = zeros(size(eta));
    En = xop2_2d(En,u.^2,0,am);
    En = yop2_2d(En,v.^2,1,am);
    En = 0.5* (g*eta.^2 + En.*depth);
    
    E(n) = dx*dy*sum(En(:));
    
    tn8(n) = n*dt; 
end

%%
Mo.E = E;
Mo.V = Vol;
Mo.t = tn8;
save('Mo8','Mo');

%%
MOALL(1).M = load('Mo2');
MOALL(2).M = load('Mo4');
MOALL(3).M = load('Mo8');

figure

subplot 211
plot(MOALL(1).M.Mova.t,MOALL(1).M.Mova.E,'b.','MarkerSize',5);
hold on 
plot(MOALL(2).M.Mova.t,MOALL(2).M.Mova.E,'r.','MarkerSize',5);
plot(MOALL(3).M.Mova.t,MOALL(3).M.Mova.E,'k.','MarkerSize',5);
hold off
xlabel 'time [ s]';
ylabel 'Energy';ylim([40 60])
grid on;
set(gca,'Fontsize',20,'FontName','Cambria');
legend('dx = 2m', 'dx = 4m', 'dx = 8m')

subplot 212
plot(MOALL(1).M.Mova.t,MOALL(1).M.Mova.V,'b.','MarkerSize',5);
hold on 
plot(MOALL(2).M.Mova.t,MOALL(2).M.Mova.V,'r.','MarkerSize',5);
plot(MOALL(3).M.Mova.t,MOALL(3).M.Mova.V,'k.','MarkerSize',5);
hold off
xlabel 'time [ s]';
ylabel 'Vol [ m^3]';ylim([196 197]);
grid on;
set(gca,'Fontsize',20,'FontName','Cambria');
legend('dx = 2m', 'dx = 4m', 'dx = 8m')

