

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
dx = 2; dy = 2;
M = floor(a/dx); % cell number at x direction
N = floor(b/dy); % cell number at y direction

xe = xmin:dx:xmax; % cell edge at x direction
ye = ymin:dy:ymax; % cell edge at y direction

xc = (xmin+dx/2):dx:(xmax - dx/2); % eta location x
yc = (ymin+dy/2):dy:(ymax - dy/2); % eta locatoin y

xev = xc; % v location x
yev = ye; % v locaiton y

%%  analysis value
cwave = sqrt(g*H);
omega = cwave *sqrt((m*pi/a)^2+(n*pi/b)^2);
Umag = alpha*g*m*pi/(omega*a);
%
t = 1:1:256;
%for n = 1:max(size(t))
for i = 1:length(xe)
    for j = 1:length(yc)
        
        ua(i,j) = Umag*sin(m*pi*xe(i)/a)*cos(n*pi*yc(j)/b)*sin(omega*t(256));
        
    end
end

figure(3)

contourf(xe,yc,ua',16);
colormap jet; shading interp;colorbar;
set(gca,'Fontsize',20,'Fontname','cambria')

%%
n=1;
Vmag = alpha*g*n*pi/omega/a;

t = 1:1:256;

for i = 1:length(xev)
    for j = 1:length(yev)
        
        va(i,j) = Vmag*cos(n*pi*xev(i)/a)*sin(n*pi*yev(j)/b)*sin(omega*t(256));
        
    end
end

figure(2)
contourf(va',10); shading interp;
zlim([-0.1 0.1]);colormap jet;
%title(['t =' num2str(n)])
pause(1/10)

%%

cwave = sqrt(g*H);
omega = cwave *sqrt((m*pi/a)^2+(n*pi/b)^2);

t = 1:1:256;

for i = 1:length(xc)
    for j = 1:length(yc)
        
        etaa(i,j) = alpha*cos(m*pi*xc(i)/a)*cos(n*pi*yc(j)/b)*cos(omega*t(256));
        
    end
end
 

figure(2)
contourf(xc,yc,etaa',10)
shading interp;
%title (['\eta t =' num2str(n)]);
colorbar; colormap jet;caxis([-0.1 0.1])
set(gca,'Fontsize',20,'Fontname','cambria')


