function [u,v,eta] = swerk3(u,v,eta, dx,dy,dt, g, depth, fcoriolis);

% stage 1
[ru, rv, rp] = swerhs(u,v,eta, dx,dy,g, depth, fcoriolis);

u1 = u + ru*dt;
v1 = v + rv*dt;
eta1 = eta + rp*dt;

% stage 2
[ru, rv, rp] = swerhs(u1,v1,eta1, dx,dy,g, depth, fcoriolis);

u2 = 0.75*u + 0.25*(u1 + dt*ru);
v2 = 0.75*v + 0.25*(v1 + dt*rv);
eta2 = 0.75*eta + 0.25*(eta1 + dt*rp);

% stage 3
[ru, rv, rp] = swerhs(u2,v2,eta2, dx,dy,g, depth, fcoriolis);

u = (u + 2.0 * (u2 + dt*ru))/3;
v = (v + 2.0 * (v2 + dt*rv))/3;
eta = (eta + 2.0 * (eta2 + dt*rp))/3;