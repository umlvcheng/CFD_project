
%% linerized shallow water equation

clear;
close all;
clc;
%%
% standing wave case

% define variable
H = 10;
a = 1280; % longth of the bay
b = 640;  % width of the bay
alpha = 0.1;
m = 2;
n = 1;
g = 10; % m/s^2;

tend = 256; % the end of intergrations

% define the domain

xmin = -a/2; xmax = a/2;
ymin = -b/2; ymax = b/2;

% define the BC of the linerized SWE
periodic = 1;

% define the compute grid

dx = 20; dy = 20; % grid spacing
%dx = 2; dy = 2;
M = a/dx; % cell number at x direction
N = b/dy; % cell number at y direction

xe = xmin:dx:xmax; % cell edge at x direction
ye = ymin:dy:ymax; % cell edge at y direction

xc = (xmin+dx/2):dx:(xmax - dx/2); % eta location x
yc = (ymin+dy/2):dy:(ymax - dy/2); % eta locatoin y

xev = xc; % v location x
yev = ye; % v locaiton y


am = [0.5  0.5];      % averaging coefficients
ad = [1/dx -1/dx];      % gradient coefficients

[fp,fu,fv,fz]=setmatrices(M,N); % define functions on C-grid points

ut = zeros(size(fu));
gz = zeros(size(fz));
vt = zeros(size(fv));
gp = zeros(size(fp));
etat = zeros(size(fp));

% Define integration time and number of time steps

Courant  = 0.5;              % Set courant number
c  = sqrt(g*H);              % unit velocity
dt = Courant*dx/c;           % Courant number is fixed
ntimesteps=2*(xmax-xmin)/c/dt; % enough steps for a SINGLE circuit

% initial condition

%u(x,y,0) = 0;
%y(x,y,0) = 0;
%eta(x,y,0) = alpha.*cos(m*pi*x/a).*cos(n*pi*y/b);

for i = 1:max(size(xc))
    for j = 1:max(size(yc))
        eta(i,j) = alpha*cos(m*pi*xc(i)/a).*cos(n*pi*yc(j)/b);
    end
end
u = zeros(size(ut));
v = zeros(size(vt));

%contourf(eta)

%%  shallow water equation

for n = 1:ntimesteps;
    
    % stage 1
    ut = -xop2_2d(ut,eta,0,ad); % tendency for stage 1
    u = u +ut*dt;  % Stage 1 update
    vt = -yop2_2d(vt,eta,0,ad);
    v = v + vt*dt;
    etat = -(xop2_2d(eta,H*u,0,ad) + yop2_2d(eta,H*v,0,ad));
    eta = eta + etat*dt;
    
    
    figure(1)
    subplot 311
    %surf(,eta);
    contour(xc,yc,eta',10);
    zlim([-0.1 0.1])
    title (['\eta t =' num2str(n)]);
    colorbar; colormap jet;caxis([-0.1 0.1])
    set(gca,'Fontsize',20)
    
    subplot 312
    %surf(u);
    contour(xe,yc,u',10);
    title 'u'
    colorbar;colormap jet;%caxis([-0.05 0.05])
    set(gca,'Fontsize',20)
    
    
    subplot 313
    %surf(v);
    contour(xev,yev,v',10);
    title 'v';colormap jet;%caxis([-0.02 0.02])
    set(gca,'Fontsize',20)
    colorbar
    
    % monitor the total volume of displace water
    eta_sum = sum(sum(eta));
    Vol(n,:) = dx*dy*eta_sum;
    
    % monitor the energy defined on the p=points
    
    uu = u.^2; vv = v.^2;
    uu = xop2_2d(eta,uu,0,am);
    vv = yop2_2d(eta,vv,0,am);
    
    Ep(n,:) = (sum(sum(g*eta.^2/2)));
    Ek(n,:) = H*(sum(sum(uu/2))) + H*(sum(sum(vv/2)));
    E(n,:) = g*(sum(sum(eta.^2/2))) + H*(sum(sum(uu/2))) + H*(sum(sum(vv/2)));
    
    pause()
    
end

figure

subplot 221
plot(E,'k.','MarkerSize',15);
set(gca,'Fontsize',20);
xlabel 'timestep';
ylabel 'Energy'

subplot 222
plot(Ep,'k.','MarkerSize',15);
set(gca,'Fontsize',20);
xlabel 'timestep';
ylabel 'Potential Energy'

subplot 223
plot(Ek,'k.','MarkerSize',15);
set(gca,'Fontsize',20);
xlabel 'timestep';
ylabel 'Kinetic Energy'


subplot 224
plot(Vol,'k.','MarkerSize',15);
set(gca,'Fontsize',20);
xlabel 'timestep';
ylabel 'Vol'


%%  analysis value
cwave = sqrt(g*H);
omega = cwave *sqrt((m*pi/a)^2+(n*pi/b)^2);
Umag = alpha*g*m*pi/omega/a;
%
t = 1:1:256;
for n = 1:max(size(t))
    for i = 1:max(size(xe))
        for j = 1:max(size(ye))
            
            ua(i,j) = Umag*sin(m*pi*xe(i)/a)*cos(n*pi*ye(j)/b)*sin(omega*t(n));
            
            
        end
    end
    
    
    %surf(ua)
    contour(ua',16);
    colormap jet
    zlim([-0.1 0.1])
    colorbar; colormap jet;
    set(gca,'Fontsize',20)
    pause(1/10)
    
end

%%
n=1;
Vmag = alpha*g*n*pi/omega/a;


t = 1:1:256;
for n = 1:max(size(t))
    for i = 1:max(size(xev))
        for j = 1:max(size(yev))
            
            va(i,j) = Vmag*cos(n*pi*xev(i)/a)*sin(n*pi*yev(j)/b)*sin(omega*t(n));
            
        end
    end
    
    figure(2)
    surf(va)
    zlim([-0.1 0.1])
    title(['t =' num2str(n)])
    pause(1/10)
    
end

%%
t = 1:1:64;
for n = 1:max(size(t))
    for i = 1:max(size(xc))
        for j = 1:max(size(yc))
            
            eta(i,j) = alpha*cos(m*pi*xc(i)/a)*cos(n*pi*yc(j)/b)*cos(omega*t(n));
            
        end
    end
    
    figure(2)
    contour(eta')
    
    title (['\eta t =' num2str(n)]);
    colorbar; colormap jet;caxis([-0.1 0.1])
    set(gca,'Fontsize',20)
    pause()
    
end



%%
for i = 1:max(size(u))
    
    
    u(i,j) = rk3(u);
    
    
end
%%
% RK3

for n = 1:256
    for i = 1:(max(size(xc))-1);
        for j = 1:(max(size(yc))-1);
            
            u(i+1,j) = -(eta(i+1,j) - eta(i,j))/dx;
            u(1,:) = 0;
            u(max(size(xc)),:) = 0;
            
            v(i,j+1) = -(eta(i,j+1) - eta(i,j))/dy;
            v(:,1) = 0;
            v(:,max(size(yc))) = 0;
            
            
            eta(i,j) = -((H*u(i+1,j) - H*u(i,j))/dx +(H*v(i,j+1) - H*v(i,j)/dy));
            eta(max(size(yc)),:) = 0;
            eta(:,max(size(yc))) = 0;
            
            
        end
    end
    
    figure(1)
    subplot 311
    surf(eta)
    
    title '\eta'
    
    subplot 312
    surf(u)
    title 'u'
    
    subplot 313
    surf(v)
    title 'v'
    
    pause()
    
end
%%



