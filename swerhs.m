function [ru, rv, rp] = swerhs(u,v,eta, dx,dy, g, depth, fcoriolis)

am =  [0.5  0.5];       % averaging coefficients
adx = [1.0 -1.0]/dx;   % x dir gradient coefficients
ady = [1.0 -1.0]/dy;   % y dir gradient coefficients

alin = 0.0;
%alin = 1.0;% for nonlinear case

ru = zeros(size(u));
rv = zeros(size(v));

thead = g*eta;
thead = xop2_2d(thead, u.^2, 1.0, alin*0.5*am);
thead = yop2_2d(thead, v.^2, 1.0, alin*0.5*am);

ru = -xop2_2d(ru,thead,0,adx); ru(1,:) = 0.0; ru(end,:) = 0.0;
rv = -yop2_2d(rv,thead,0,ady); rv(:,1) = 0.0; rv(:,end) = 0.0;

%h = depth + alin * eta;
h = depth;

hbx = zeros(size(u));
hbx = xop2_2d(hbx,h, 0, am);hbx(1,:) = h(1,:);hbx(end,:) = h(end,:);
U = hbx.*u;

hby = zeros(size(v));
hby = yop2_2d(hby,h, 0, am); hby(:,1) = h(:,1); hby(:,end) = h(:,end);
V = hby.*v;

hbxy = zeros(size(fcoriolis));
hbxy = xop2_2d(hbxy, hby, 0.0, am); hbxy(1,:) = hby(1,:); hbxy(end,:) = hby(end,:);

q = fcoriolis;
q = xop2_2d(q, v, 1.0, alin*adx);
q = yop2_2d(q, u, 1.0,-alin*ady);
q = q./hbxy;

EnstrophyBudget = 0.5*dx*dy* sum( sum(hbxy .* q.^2) );

rp = -(xop2_2d(eta,U,0,adx) + yop2_2d(eta,V,0,ady));


