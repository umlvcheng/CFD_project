
%% linerized shallow water equation

clear;
close all;
clc;
% standing wave case

% define variable
%H = 10;
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

dx = 4; dy = 4; % grid spacing
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

depth = ones([M N]);
Hu = 10 + 5 * tanh((xc - 100)/16);
for i = 1 : max(size(yc))
    depth(:,i) = Hu;
end

Courant  = 0.5;                % Set courant number
c  = sqrt( g*max( depth(:)));  % maximum gravity wave speed
dt = Courant*dx/c;             % Courant number is fixed
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

for n = 1:ntimesteps+1;
    
    [u,v,eta] = swerk3(u,v,eta, dx,dy,dt, g, depth, fcoriolis);
    
%     if floor(n * dt) == 32 ;
%         sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
%         pause();
%     elseif floor(n * dt) == 64 ;
%         sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
%         pause();
%     elseif  floor(n * dt) == 128 ;
%         sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
%         pause();
%     else
%     if  floor(n * dt) == 256 ;
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
Mova.E = E;
Mova.V = Vol;
Mova.t = tn8;
save('Mo4','Mova');

%%
figure

subplot 211
plot(E,'k.','MarkerSize',15);
xlabel 'timestep';
ylabel 'Energy';ylim([0 100])
grid on;
set(gca,'Fontsize',20,'FontName','Cambria');

subplot 212
plot(Vol,'k.','MarkerSize',15);
xlabel 'timestep';
ylabel 'Vol [ m^3]';ylim([100 300]);
grid on;
set(gca,'Fontsize',20,'FontName','Cambria');

