% Programming example for finite differences in matlab
% Solve the problem p_t + u p_x = 0
%

% Define the domain:
xmin =-1; xmax=1;

% Define the BC of the problem
periodic=1;

% Define the computational grid;
N=64;              % number of intervals
dx=(xmax-xmin)/N;  % grid spacing
xe=xmin:dx:xmax;   % cell edges

% Define integration time and number of time steps
Courant  = 0.5;              % Set courant number
u  = 1.0;                    % unit velocity
dt = Courant*dx/u;           % Courant number is fixed
ntimesteps=(xmax-xmin)/u/dt; % enough steps for a SINGLE circuit

% define initial condition function
%init = inline('exp(-sin(pi*x).^2)');
 alpha = 4;
 init = inline('tanh(4*(cos(pi*x) + 1/2)) + exp(-(cos(pi*x).^2))');
%init = inline('max(0,1-4*x.^2)');

% Initial Conditions
p = init(xe);

umin = min(p); umax=max(p); du = (umax-umin)*0.1;
axislimits=[xmin xmax min(p)-du max(p)+du];
plot(xe,p,'k');
axis(axislimits);
disp('Initial Conditions. Hit return to continue');
pause

pexact = p; % save initial conditions for later use
isnap=4;     % Define how often (# of time steps) to plot solution

for it = 1:ntimesteps

  p = rk3(p,u,dx,dt,periodic); % RK3 time step

  time = it*dt;
%%%%%%%%%%%%%Plot numerical solution %%%%%%%%%%%%%%%%%%%%%
  if (mod(it,isnap)==0)
    plot(xe,p,'k');
    legend('Numerical');
    xlabel('x'); ylabel('u');
    axis(axislimits);

    disp('Hit return to continue');
    
    pause
  end
%%%%%%%%%%%%%Plot numerical solution %%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%Plot final solution %%%%%%%%%%%%%%%%%%%%%
plot(xe,pexact,'r', xe,p,'k','LineWidth',3);
legend('Exact','Numerical');
xlabel('x'); ylabel('p');
axis(axislimits);
set(gca,'FontName','Cambira','FontSize',16)
%%%%%%%%%%%%%Error Measures%%%%%%%%%%%%%%%%%%%%%
err     = (p-pexact);
emax=max(abs(err));
erms=sqrt(sum(err.^2))*dx;
disp(strcat('Errors:', num2str([N emax erms]) ) );

%%
text(-0.8,-0.1,'N = 64, Emax = 0.117894, Erms = 0.01078','Fontsize',16)
text(-0.8,-0.1,'N = 128, Emax = 0.06552023, Erms = 0.003680782','Fontsize',16)

text(-0.8,-0.1,'N = 128, Emax = 0.1638743, Erms = 0.009361916','Fontsize',16)
text(-0.8,-0.1,'N = 64, Emax = 0.37711, Erms = 0.0373239','Fontsize',16)

%%

