
%% linerized shallow water equation

clear;
close all;
clc;
% standing wave case

% define variable
H = 10;
a = 1280; % longth of the bay
b = 640;  % width of the bay
alpha = 0.1;m = 2;n = 1;g = 10; % m/s^2;

% define the domain
xmin = -a/2; xmax = a/2;
ymin = -b/2; ymax = b/2;

% define the BC of the linerized SWE
periodic = 1;

% define the compute grid
%dx = 10; dy = 10; % grid spacing
dx = 20; dy = 20;
M = floor(a/dx); % cell number at x direction
N = floor(b/dy); % cell number at y direction

xe = xmin:dx:xmax; % cell edge at x direction
ye = ymin:dy:ymax; % cell edge at y direction

xc = (xmin+dx/2):dx:(xmax - dx/2); % eta location x
yc = (ymin+dy/2):dy:(ymax - dy/2); % eta locatoin y

xev = xc; % v location x
yev = ye; % v locaiton y


am = [0.5  0.5];      % averaging coefficients
ad = [1/dx -1/dx];      % gradient coefficients

% Define integration time and number of time steps
depth = H*ones([M N]);
Courant  = 0.5;              % Set courant number
c  = sqrt( g*max( depth(:) ) );              % maximum gravity wave speed
dt = Courant*dx/c;           % Courant number is fixed
ntimesteps=2*(xmax-xmin)/c/dt; % enough steps for a SINGLE circuit

% initial condition

u = zeros([M+1 N]);
v = zeros([M N+1]);
eta = zeros([M N]);

for j = 1:length(yc)
    eta(:,j) = alpha*cos(m*pi*xc(:)/a).*cos(n*pi*yc(j)/b);
end

fcoriolis = zeros([M+1 N+1]);
n=0;
sweplots(u,v,eta,xc,yc,xe,ye,n,dt,dx,dy);

%% shallow water equation

Vol = zeros([ntimesteps 1]);
E = zeros([ntimesteps 1]);

for n = 1:ntimesteps;
    
    [u,v,eta] = swerk3(u,v,eta, dx,dy,dt, g, depth, fcoriolis);

    if n * dt == 32 ;
        sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
        pause();
    elseif n * dt == 64 ;
        sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
        pause();
    elseif n * dt == 128 ;
        sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
        pause();
    elseif n * dt == 256 ;
        sweplots(u,v,eta,xc,yc,xe,ye,n,dt);
        pause();
        
    end
    % monitor the total volume of displace water
    Vol(n) = dx*dy*sum(eta(:));
    
    % monitor the energy defined on the p=points
    En = zeros(size(eta));
    En = xop2_2d(En,u.^2,0,am);
    En = yop2_2d(En,v.^2,1,am);
    En = 0.5* (g*eta.^2 + En.*H);
    
    E(n) = dx*dy*sum(En(:));
    tn20(n) = n*dt;
    
end

%%
Mostd.E = E;
Mostd.V = Vol;
Mostd.t = tn20;
save('Mo20','Mostd');

%%
MOALL(1).M = load('Mo2');
MOALL(2).M = load('Mo5');
MOALL(3).M = load('Mo10');
MOALL(4).M = load('Mo20');

figure

subplot 211
plot(MOALL(1).M.Mostd.t,MOALL(1).M.Mostd.E,'b.','MarkerSize',5);
hold on 
plot(MOALL(2).M.Mostd.t,MOALL(2).M.Mostd.E,'r.','MarkerSize',5);
plot(MOALL(3).M.Mostd.t,MOALL(3).M.Mostd.E,'k.','MarkerSize',5);
plot(MOALL(4).M.Mostd.t,MOALL(4).M.Mostd.E,'m.','MarkerSize',5);
hold off
xlabel 'time [ s]';
ylabel 'Energy';ylim([1.022*1e4 1.025*1e4])
grid on;
set(gca,'Fontsize',20,'FontName','Cambria');
legend('dx = 2m', 'dx = 5m', 'dx = 10m','dx = 20m')

subplot 212
plot(MOALL(1).M.Mostd.t,MOALL(1).M.Mostd.V,'b.','MarkerSize',5);
hold on 
plot(MOALL(2).M.Mostd.t,MOALL(2).M.Mostd.V,'r.','MarkerSize',5);
plot(MOALL(3).M.Mostd.t,MOALL(3).M.Mostd.V,'k.','MarkerSize',5);
plot(MOALL(4).M.Mostd.t,MOALL(4).M.Mostd.V,'m.','MarkerSize',5);
hold off
xlabel 'time [ s]';
ylabel 'Vol [ m^3]';ylim([-1e-12 5e-12]);
grid on;
set(gca,'Fontsize',20,'FontName','Cambria');
legend('dx = 2m', 'dx = 5m', 'dx = 10m','dx = 20m')
%% monitoring the rms errors of the velocity and pressure variables at time t = 256;

rms_eta = sum(sum( eta - etaa))/M/N;
rms_u = sum(sum(u - ua))/(M+1)/N;
rms_v = sum(sum(v -va))/M/(N+1);





