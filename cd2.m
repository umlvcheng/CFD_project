function du = cd2(ue,dx,periodic)
%
% Returns the 2nd-order centered difference derivative of the function 
% defined in the array ue
% Usage is
% du = cd2(ue,dx,periodic)
% where 
%  dx        is the grid spacing (assumed uniform)
%  ue        is the array to be differentiated 
%  periodic  is a flag to indicate whether to use periodicity to deal with boundary ghost points
%  du        is the array containing the finite difference approximation
%
  Ne = length(ue); % number of edges
  du(2:Ne-1) = (ue(3:Ne)-ue(1:Ne-2))/(2*dx);     % Centered 2nd order
  if (periodic==1)
    du(1) = (ue(2)-ue(Ne-1))/(2*dx);     % Centered 2nd order
    du(Ne)= (ue(2)-ue(Ne-1))/(2*dx);
  end
