function f = init(x);
%
% Initializes a function
%
  f = ones(size(x));
  ind = find(abs(x)>0.5);
  f(ind)=0.0;
